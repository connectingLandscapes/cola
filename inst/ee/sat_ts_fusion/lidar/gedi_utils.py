import os
import ee
import ee.cli.commands
import requests
import zipfile
import geopandas as gpd
import pandas as pd
import shapely.geometry as sg
from google.cloud import storage
import subprocess
from pathlib import Path
import numpy as np
from pandas import DataFrame


def create_geo_grid(xmin:float = -180, xmax:float = 180,
                    ymin:float = -90, ymax:float = 90,
                    res:float = 1):
    """Create a regularly spaced polygon grid with geographic coordinates

    :param xmin: minimum longitude in decimal degress
    :param xmax: maximum longitude in decimal degress
    :param ymin: minimum latitude in decimal degress
    :param ymax: maximum latitude in decimal degress
    :param res: resolution of the grid in decimal degrees
    :return: geopandas data frame
    """

    grid_cells = []
    for lon in range(xmin, xmax, res):
        for lat in range(ymin, ymax, res):
            # Create a polygon representing the grid cell
            cell = sg.Polygon([
                (lon, lat),
                (lon + res, lat),
                (lon + res, lat + res),
                (lon, lat + res)
            ])
            grid_cells.append(cell)

    # Convert to GeoDataFrame
    grid_gdf = gpd.GeoDataFrame(geometry=grid_cells, crs="EPSG:4326")

    return grid_gdf


def create_grid_covering_aoi(aoi_shape: str = None,
                             crs: str = None,
                             width: float = 10000):
    """Create a grid of square polygons covering an AOI shapefile

    :param aoi_shape: path and name of the AOI shapefile
    :param crs: coordinate reference system specified as EPSG code, like EPSG:3310
    :param width: the square polygon width in units of the crs
    :return: a geo dataframe of square polygon grids
    """

    aoi_gdf = gpd.read_file(aoi_shape)
    aoi_gdf_rpj = aoi_gdf.to_crs(crs=crs)
    xmin_rpj, ymin_rpj, xmax_rpj, ymax_rpj = aoi_gdf_rpj.total_bounds

    cols = list(np.arange(xmin_rpj, xmax_rpj + width, width))
    rows = list(np.arange(ymin_rpj, ymax_rpj + width, width))

    polygons = []
    for x in cols[:-1]:
        for y in rows[:-1]:
            polygons.append(sg.Polygon([(x, y), (x + width, y), (x + width, y + width), (x, y + width)]))

    grid = gpd.GeoDataFrame({'geometry': polygons}).set_crs(crs)

    # select intersecting grids and clip to the AOI
    selected_grids = grid[grid.geometry.intersects(aoi_gdf_rpj.union_all())]
    clipped_grids = gpd.clip(selected_grids, aoi_gdf_rpj)
    # assign a grid id
    clipped_grids["grid_id"] = clipped_grids.index

    return clipped_grids


def decimate_with_grid(gdf_grid: gpd.GeoDataFrame = None,
                       gdf_points: gpd.GeoDataFrame = None,
                       percentile_threshold: float = 0.2,
                       random_seed: int = 123):
    """Decimate (reduce) the number of GEDI shots in a regular grid to create a more even spatial sample

    :param gdf_grid: grid polygon geo dataframe
    :param gdf_points: points geo dataframe
    :param percentile_threshold: the point density quantile threshold to use for subsampling
    :param random_seed: for reproducibility, defaults to 123
    :return: a geo dataframe of GEDI shots that has been subsampled per grid
    """

    # add area column, ideally this is projected area
    gdf_grid["area"] = gdf_grid.area

    # spatially join the grid and points
    joined_gdf = gpd.sjoin(gdf_points, gdf_grid.to_crs('EPSG:4326'), how='left', predicate='within')

    # estimate shot density per grid
    shots_per_grid = joined_gdf.groupby("grid_id").size().rename("shot_count")
    gdf_grid["shot_count"] = shots_per_grid
    gdf_grid["shot_count"] = gdf_grid["shot_count"].fillna(0).astype(int)
    gdf_grid["density"] = gdf_grid["shot_count"] / gdf_grid["area"]

    # compute shot density threshold to use for subsampling
    gdf_grid_f = gdf_grid[gdf_grid["shot_count"] >= 20]
    density_thresh = gdf_grid_f["density"].quantile(percentile_threshold)
    gdf_grid["shot_sample_n"] = (density_thresh * gdf_grid["area"]).round()
    # for the grids that fall below the density threshold, make sure not to sample more shots than are available
    gdf_grid.loc[gdf_grid['shot_sample_n'] >= gdf_grid['shot_count'], 'shot_sample_n'] = gdf_grid['shot_count']

    # subsample based on the density threshold
    gdf_points_ss_list = []
    for grid_id, group in joined_gdf.groupby('grid_id'):
        samp_size = gdf_grid.loc[gdf_grid['grid_id'] == grid_id, 'shot_sample_n'].iloc[0]

        # Ensure we don't try to sample more points than are available in the group
        act_samp_size = int(min(samp_size, len(group)))

        if act_samp_size > 0:
            gdf_points_ss_grid = group.sample(n=act_samp_size,
                                              random_state=random_seed)  # Use random_state for reproducibility
            gdf_points_ss_list.append(gdf_points_ss_grid)

    gdf_points_ss_merge = gpd.GeoDataFrame(pd.concat(gdf_points_ss_list, ignore_index=True), crs=gdf_points.crs)

    # Drop the columns
    columns_to_drop = ['geometry', 'area', 'index_right']
    gdf_points_ss_merge_drop = gdf_points_ss_merge.drop(columns=columns_to_drop)

    return gdf_points_ss_merge_drop


def list_l24a_zips(base_url:str = 'https://rcdata.nau.edu/geode_data/GEDIv002_L0204A_20190417to20230316_proc202312/tables/',
                   aoi_shape:str = None):
    """List GEDI L2L4A zipped files that overlap with an area of interest

    :param base_url: the URL where the zipped tables are hosted
    :param aoi_shape: shapefile covering the area of interest
    :return: a list of zipped table URLs
    """

    # intersect aoi with 1 degree grid
    grid_1deg_gdf = create_geo_grid(xmin=-180, xmax=180, ymin=-52, ymax=52, res=1)
    aoi_gdf = gpd.read_file(aoi_shape)
    aoi_gdf_rpj = aoi_gdf.to_crs(crs="EPSG:4326")
    selected_grids = grid_1deg_gdf[grid_1deg_gdf.geometry.intersects(aoi_gdf_rpj.unary_union)]
    def get_tile_url(row):
        xmin = int(row.geometry.bounds[0])
        ymin = int(row.geometry.bounds[1])

        if (xmin < 0):
            ew = "w"
        else:
            ew = "e"

        if (ymin < 0):
            ns = "s"
        else:
            ns = "n"

        # Construct the tile ID
        ew_pad = f"{abs(xmin):03d}"+ew
        ns_pad = f"{abs(ymin):02d}"+ns
        tileid= "g"+ew_pad+ns_pad

        url = base_url+"gediv002_l2l4a_va_"+tileid+".zip"
        return url

    selected_grids["url"] = selected_grids.apply(get_tile_url, axis=1)
    selected_grids_list = selected_grids["url"].tolist()
    list_len = len(selected_grids_list)
    print(f"{list_len} GEDI quality-filtered tables will be downloaded from {base_url}...\n")

    return selected_grids["url"].tolist()


def clip_l2l4a_table(table_df:pd.DataFrame = None,
                     aoi_shape:str = None):
    """Clip GEDI L2L4A table

    :param table_df: pandas DataFrame with x,y geographic coordinates
    :param aoi_shape: the shapefile to use for clipping the table
    :return: a table clipped to a shapefile
    """

    aoi_gdf = gpd.read_file(aoi_shape)
    aoi_gdf_rpj = aoi_gdf.to_crs(crs="EPSG:4326")
    table_df["geometry"] = [sg.Point(x, y) for x, y in zip(table_df["lon_lm_a0"], table_df["lat_lm_a0"])]
    gdf = gpd.GeoDataFrame(table_df, geometry="geometry", crs="EPSG:4326")
    gdf_clip = gdf.clip(aoi_gdf_rpj)

    return gdf_clip


def filter_l2l4a_by_years_doys(df:pd.DataFrame = None,
                               start_year: int = 2019,
                               end_year: int = 2023,
                               start_doy: int = 1,
                               end_doy: int = 366):
    """

    :param df: GEDI L2L4A table as a pandas DataFrame
    :param start_year: first year to include
    :param end_year: last year to include
    :param start_doy: first day of year to include
    :param end_doy: last day of year to include
    :return: a filtered pandas DataFrame
    """

    df_f = df[(df['date_dec'] >= start_year) & (df['date_dec'] <= end_year) & (df['doy'] >= start_doy) & (df['doy'] <= end_doy)]

    return df_f



def filter_l2l4a_by_elevation_monthly(df:pd.DataFrame = None,
                                      max_elev_ls:list = None):
    """Filter a DataFrame of GEDI shots by elevation on a monthly basis

    :param df: a pandas DataFrame
    :param max_elev_ls: a list of the maximum allowable elevations of GEDI shots per month (January to December)
    :return: a filtered pandas DataFrame
    """
    df_jan = df[(df['elev_lm_a0'] <= max_elev_ls[0]) & (df['doy'] >= 1) & (df['doy'] <= 31)]
    df_feb = df[(df['elev_lm_a0'] <= max_elev_ls[1]) & (df['doy'] >= 32) & (df['doy'] <= 59)]
    df_mar = df[(df['elev_lm_a0'] <= max_elev_ls[2]) & (df['doy'] >= 60) & (df['doy'] <= 90)]
    df_apr = df[(df['elev_lm_a0'] <= max_elev_ls[3]) & (df['doy'] >= 91) & (df['doy'] <= 120)]
    df_may = df[(df['elev_lm_a0'] <= max_elev_ls[4]) & (df['doy'] >= 121) & (df['doy'] <= 151)]
    df_jun = df[(df['elev_lm_a0'] <= max_elev_ls[5]) & (df['doy'] >= 152) & (df['doy'] <= 181)]
    df_jul = df[(df['elev_lm_a0'] <= max_elev_ls[6]) & (df['doy'] >= 182) & (df['doy'] <= 212)]
    df_aug = df[(df['elev_lm_a0'] <= max_elev_ls[7]) & (df['doy'] >= 213) & (df['doy'] <= 243)]
    df_sep = df[(df['elev_lm_a0'] <= max_elev_ls[8]) & (df['doy'] >= 244) & (df['doy'] <= 273)]
    df_oct = df[(df['elev_lm_a0'] <= max_elev_ls[9]) & (df['doy'] >= 274) & (df['doy'] <= 304)]
    df_nov = df[(df['elev_lm_a0'] <= max_elev_ls[10]) & (df['doy'] >= 305) & (df['doy'] <= 334)]
    df_dec = df[(df['elev_lm_a0'] <= max_elev_ls[11]) & (df['doy'] >= 335) & (df['doy'] <= 366)]

    df_filt = pd.concat([df_jan, df_feb, df_mar, df_apr, df_may, df_jun, df_jul, df_aug, df_sep, df_oct, df_nov, df_dec], axis=0)

    return df_filt


def dl_unzip_filter_clip_merge_l24a_zips(zip_url_list:list = None,
                                         dl_dir:str = None,
                                         start_year: int = 2019,
                                         end_year: int = 2023,
                                         start_doy: int = 1,
                                         end_doy: int = 366,
                                         max_elev_ls: list = None,
                                         aoi_shape:str = None,
                                         gedi_keep_cols:list = None,
                                         max_shots_per_tile: int = None,
                                         decimate:bool = True,
                                         grid_crs:str = None,
                                         grid_width: int = 10000,
                                         output_table_path:str = None):
    """Download and unzip GEDI L2L4A zipped tables

    :param zip_url_list: list of zip files to be downloaded
    :param dl_dir: the directory to store the downloaded zip files
    :param start_year
    :param end_year
    :param start_doy
    :param end_doy
    :param max_elev_ls
    :param aoi_shape: the shapefile to use for clipping the table
    :param gedi_keep_cols:
    :param max_shots_per_tile
    :param grid_crs
    :param grid_width
    :param output_table_path: the path and file name for the merged csv table
    :return: unzipped tables clipped to a shapefile and saved as CSVs
    """

    dl_zip_dir = dl_dir+'0_zips/'
    tables_dir = dl_dir+'1_tables/'
    os.makedirs(dl_zip_dir, exist_ok=True)
    os.makedirs(tables_dir, exist_ok=True)

    # download, extract, and clip all zip files
    # save resulting geopandas df in a list
    gdf_list = []
    i = 1
    list_len = len(zip_url_list)
    for url in zip_url_list:
        zip_filename = os.path.join(dl_zip_dir, os.path.basename(url))
        path_zip = Path(zip_filename)
        table_base_noext = path_zip.stem
        print(f"Working on table {i} out of {list_len} - {table_base_noext}...")

        # download
        if os.path.exists(zip_filename):
            print(f"...downloaded zip file exists")
        else:
            response = requests.get(url)
            with open(zip_filename, "wb") as f:
                print("...downloading zip file")
                f.write(response.content)

        # unzip
        table_filename = tables_dir + table_base_noext + '.csv'
        if os.path.exists(table_filename):
            print(f"...unzipped table exists")
        else:
            with zipfile.ZipFile(zip_filename, "r") as zip_ref:
                print("...unzipping table")
                zip_ref.extractall(tables_dir)

        # read as pandas dataframe
        print("...reading table as geo dataframe")
        df = pd.read_csv(table_filename, usecols=gedi_keep_cols)

        # filter temporally
        df_f = filter_l2l4a_by_years_doys(df=df, start_year=start_year, end_year=end_year, start_doy=start_doy, end_doy=end_doy)

        # filter by elevation on a monthly basis
        df_f = filter_l2l4a_by_elevation_monthly(df=df_f, max_elev_ls=max_elev_ls)

        # clip to the region of interest
        if aoi_shape is not None:
            print("...clipping geo dataframe to AOI")
            gdf = clip_l2l4a_table(table_df=df_f, aoi_shape=aoi_shape)
            #print(f"...saving {clipped_table_filename} \n")
            #clipped_table.to_csv(clipped_table_filename, index=False)
        else:
            df_f["geometry"] = [sg.Point(x, y) for x, y in zip(df_f["lon_lm_a0"], df["lat_lm_a0"])]
            gdf = gpd.GeoDataFrame(df, geometry="geometry", crs="EPSG:4326")

        if max_shots_per_tile is not None:
            gdf_len = len(gdf)
            if gdf_len > max_shots_per_tile:
                print(f"...randomly sampling {max_shots_per_tile} rows from geo dataframe with {gdf_len} rows")
                gdf_samp = gdf.sample(max_shots_per_tile)
                gdf_list.append(gdf_samp)
            else:
                gdf_list.append(gdf)

        else:
            gdf_list.append(gdf)

        os.remove(table_filename)
        i = i+1

    # merge the list of geo dataframes
    print(f"\nMerging geodataframes...")
    gdf_merged = gpd.GeoDataFrame(pd.concat(gdf_list, ignore_index=True))
    gdf_merged['random_1'] = np.random.rand(len(gdf_merged))
    gdf_merged_len = len(gdf_merged)
    print(f"...{gdf_merged_len} rows")

    if decimate:
        print(f"\nDecimating geodataframe...")
        # make a covering grid
        grid = create_grid_covering_aoi(aoi_shape=aoi_shape, crs=grid_crs, width = grid_width)
        # decimate shots using the grid
        gdf_merged_out = decimate_with_grid(gdf_grid=grid, gdf_points=gdf_merged)
        gdf_merged_out_len = len(gdf_merged_out)
        print(f"...{gdf_merged_out_len} rows")
    else:
        gdf_merged_out = gdf_merged

    if output_table_path is not None:
        print(f"\nSaving table - {output_table_path} \n")
        gdf_merged_out.to_csv(output_table_path, index=False)



def upload_file_to_gcs(local_path_to_file:str = None,
                       gcs_bucket_name:str = None,
                       gcs_dest_file:str = None,
                       service_account_key:str = None):
    """Upload a file stored locally to a Google Cloud Storage bucket

    :param local_path_to_file: the path + file name of the local file to be uploaded
    :param gcs_bucket_name: the bucket name
    :param gcs_dest_file: the path+file name in the destination bucket (do NOT include gs://)
    :param service_account_key: the path + file name of the local JSON service account key file
    :return: Google Cloud Storage path + file name
    """
    # Explicitly use service account credentials by specifying the private key
    # file.
    storage_client = storage.Client.from_service_account_json(service_account_key)

    bucket = storage_client.get_bucket(gcs_bucket_name)
    blob = bucket.blob(gcs_dest_file)
    gcs_path = 'gs://' + gcs_bucket_name + '/' + gcs_dest_file
    print(f"Uploading merged table to GCS - {gcs_path}\n")
    blob.upload_from_filename(local_path_to_file)
    return gcs_path


def upload_point_csv_to_gee_fc(gcs_csv_path:str = None,
                               fc_asset_id:str = None,
                               gee_project_name:str = None,
                               x_col:str = "x",
                               y_col:str = "y",
                               crs:str = "EPSG:4326"):
    """Upload a csv file stored in a Google Cloud Storage bucket to GEE as a FeatureCollection

    :param gcs_csv_path: Google Cloud Storage CSV file path (should start with gs://)
    :param fc_asset_id: asset ID to use for the GEE FeatureCollection
    :param x_col: column name with x coordinates
    :param y_col: column name with y coordinates
    :param crs: The default CRS code or WKT string specifying the coordinate reference system of any geometry without one
    :return:
    """
    print(f"Uploading merged table to GEE asset- {fc_asset_id}\n")
    command1 = "earthengine set_project "+ gee_project_name
    subprocess.run(command1)
    command2 = " ".join(["earthengine upload table", '--crs', crs, '--x_column', x_col, '--y_column', y_col, '--asset_id', fc_asset_id, gcs_csv_path])
    #print(command2)
    result = subprocess.run(command2, capture_output=True, text=True)

    return fc_asset_id


def make_gedi_l2l4a_fc(start_year: int = 2019,
                       end_year: int = 2023,
                       start_doy: int = 1,
                       end_doy: int = 366,
                       max_elev_ls: list = None,
                       aoi_shape:str = None,
                       gedi_keep_cols:list = None,
                       max_shots_per_tile:int = None,
                       grid_crs:str = None,
                       grid_width: int = 10000,
                       local_working_dir:str = None,
                       gcs_bucket_name:str = None,
                       gcs_dir:str = None,
                       service_account_key:str = None,
                       gee_project_name:str = None,
                       gee_fc_asset_dir:str = None,
                       gee_fc_asset_name:str = None):
    """Download quality-filtered GEDI tables, clip to a region, decimate, and upload to GEE as a ee.FeatureCollection

    :param start_year: first year of GEDI data to use
    :param end_year: last year of GEDI data to use
    :param start_doy: first day of acceptable day of year window
    :param end_doy: last day of acceptable day of year window
    :param max_elev_ls: a list of the maximum elevation to accept on a monthly basis
    :param aoi_shape: shapefile of the area of interest
    :param gedi_keep_cols: list of GEDI column names to keep
    :param max_shots_per_tile: maximum number of shots to keep per 1x1 degree tile
    :param grid_crs: the coordinate reference system to use for defining a regularly spaced grid
    :param grid_width: the grid cell width
    :param local_working_dir: where to download and process files locally
    :param gcs_bucket_name: Google Cloud Storage bucket name
    :param gcs_dir: Google Cloud Storage bucket subdirectory
    :param service_account_key: Google Cloud project service account key file path
    :param gee_project_name: Google Cloud project name
    :param gee_fc_asset_dir: Google Earth Engine asset path
    :param gee_fc_asset_name: Google Earth Engine asset name
    :return: an ee.FeatureCollection of GEDI shots for a region
    """
    # find zipped tables that intersect AOI
    url_list = list_l24a_zips(base_url='https://rcdata.nau.edu/geode_data/GEDIv002_L0204A_20190417to20230316_proc202312/tables/',
                              aoi_shape=aoi_shape)

    # download, filter, clip and merge tables
    merged_dir = local_working_dir + '2_merged_table/'
    os.makedirs(merged_dir, exist_ok=True)
    merged_csv_file = gee_fc_asset_name + '.csv'
    dl_unzip_filter_clip_merge_l24a_zips(zip_url_list=url_list,
                                         dl_dir=local_working_dir,
                                         start_year=start_year,
                                         end_year=end_year,
                                         start_doy=start_doy,
                                         end_doy=end_doy,
                                         max_elev_ls=max_elev_ls,
                                         aoi_shape=aoi_shape,
                                         gedi_keep_cols=gedi_keep_cols,
                                         max_shots_per_tile=max_shots_per_tile,
                                         grid_crs=grid_crs,
                                         grid_width=grid_width,
                                         output_table_path=merged_dir + merged_csv_file)

    # upload the merged CSV to Google Cloud Storage bucket
    gcs_csv_file = gcs_dir + merged_csv_file
    gcs_csv_path = upload_file_to_gcs(gcs_dest_file=gcs_csv_file,
                                      local_path_to_file=merged_dir+merged_csv_file,
                                      gcs_bucket_name=gcs_bucket_name,
                                      service_account_key=service_account_key)

    # upload the merged CSV to GEE as a FeatureCollection
    gee_fc_asset_id = gee_fc_asset_dir + gee_fc_asset_name
    upload_point_csv_to_gee_fc(gcs_csv_path=gcs_csv_path, x_col="lon_lm_a0", y_col="lat_lm_a0", crs="EPSG:4326",
                               fc_asset_id=gee_fc_asset_id, gee_project_name=gee_project_name)

    return ee.FeatureCollection(gee_fc_asset_id)


def merge_l24a_tables(tables_dir:str = None,
                      merged_table_file:str = None):
    """Clip GEDI L2L4A tables to an area of interest shapefile and merge into one table

    :param tables_dir: directory with unzipped GEDI tables
    :param merged_table_file: file path for saving the merged table (optional)
    :return: merged data frame, optionally saved as CSV
    """
    print(f"Merging tables in {tables_dir}...")
    # read spatial tables and combine
    df_list = []

    i = 1
    tables_list = os.listdir(tables_dir)
    tables_list_len = len(tables_list)
    for file in tables_list:
        print(f"Working on merging table {i} out of {tables_list_len}...")
        file_path = os.path.join(tables_dir, file)
        df = pd.read_csv(file_path)
        df_list.append(df)
        i = i + 1

    merged_df = pd.concat(df_list)

    if merged_table_file is not None:
        print(f"...saving {merged_table_file} \n")
        merged_df.to_csv(merged_table_file, index=False)

    return merged_df


def ee_fc_decimate_points_wgrid(fc: ee.FeatureCollection = None,
                                region_geom: ee.Geometry = None,
                                grid_res: float = 0.1,
                                grid_crs: str = 'EPSG:4326',
                                grid_limit: int = 100,
                                export_fc_asset_id: str = None):

    # make a covering grid
    grid_fc = region_geom.coveringGrid(proj=grid_crs, scale = grid_res)

    # add a random ID column
    fc_r = fc.randomColumn(columnName='rand_id')

    def limit_pts_in_grid(grid):
        return ee.FeatureCollection(fc_r.filterBounds(grid.geometry()).limit(maximum=grid_limit, prop='rand_id', ascending=True))

    dec_fc = grid_fc.map(limit_pts_in_grid).flatten()

    if export_fc_asset_id is not None:
        task = ee.batch.Export.table.toAsset(collection=dec_fc,
                                             description='FC decimation toAsset -> '+export_fc_asset_id,
                                             assetId=export_fc_asset_id)
        task.start()

    return dec_fc

