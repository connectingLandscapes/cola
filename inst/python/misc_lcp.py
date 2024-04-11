# -*- coding: utf-8 -*-
"""
Created on Tue Jul  4 13:29:55 2023

@author: pj276
"""
import osgeo
import rasterio as rio
from rasterio.crs import CRS

#asciiToGeoTiff(r"C:\Users\pj276\Projects\UNICOR_arch\unicor\resist_sabah_example_XY200.csv.kdepaths",r"C:\Users\pj276\Projects\UNICOR_arch\unicor\resist_sabah_example_XY200.tif",crs="ESRI:102028")
asciiToGeoTiff(r"C:\Users\pj276\Projects\UNICOR_arch\unicor\resist_sabah_example_XY200.csv.addedpaths.txt",r"C:\Users\pj276\Projects\UNICOR_arch\unicor\sabah_example_XY200_unicor_crk.tif",crs="ESRI:102028")
asciiToGeoTiff(r"C:\Users\pj276\Projects\UNICOR_arch\unicor\size6_side1000px_1000000totpix.rsg", r"C:\Users\pj276\Projects\UNICOR_arch\unicor\size6_side1000px_1000000totpix.tif", crs="ESRI:102028")


cf.asciiToGeoTiff(r"C:\Users\pj276\Projects\UNICOR_arch\unicor\size6_side1000px_1000000totpix_sabah_100.csv.kdepaths",r"C:\Users\pj276\Projects\UNICOR_arch\unicor\sabah_example_100_unicor_crk_kvol.tif",crs="ESRI:102028")
cf.asciiToGeoTiff(r"C:\Users\pj276\Projects\UNICOR_arch\unicor\size6_side1000px_1000000totpix_sabah_100.csv.kdepaths",r"C:\Users\pj276\Projects\UNICOR_arch\unicor\sabah_example_100_unicor_crk_kvol2.tif",crs="ESRI:102028")

cf.asciiToGeoTiff(r"C:\Users\pj276\Projects\UNICOR_arch\unicor\resist_sabah_example_sabah_100.csv.kdepaths",r"C:\Users\pj276\Projects\UNICOR_arch\unicor\resist_sabah_example_sabah_100_unicrk.tif",crs="ESRI:102028")

cf.asciiToGeoTiff(r"C:\Users\pj276\Projects\UNICOR_arch\unicor\resist_sabah_example_sabah_2pts.kdepaths",r"C:\Users\pj276\Projects\UNICOR_arch\unicor\resist_sabah_example_sabah_2_unicrk.tif",crs="ESRI:102028")

cf.asciiToGeoTiff(r"C:\Users\pj276\Projects\UNICOR_arch\unicor\resist_sabah_example_sabah_4pts.kdepaths",r"C:\Users\pj276\Projects\UNICOR_arch\unicor\resist_sabah_example_sabah_4_unicrk.tif",crs="ESRI:102028")


# Raster name
rname = "C:/Users/pj276/Projects/UNICOR_arch/unicor/resist_sabah_example.rsg"

# Set raster projection and save to file
# Read raster and specify projection (it's not specified in the rsg formatted file)
with rio.open(rname, 'r+') as src:
    src.crs = CRS.from_string("ESRI:102028")
    profile = src.profile
    data = src.read()

# Update driver to geotiff
profile.update(driver="GTiff")
# Set output path
pathout = "C:/Users/pj276/Projects/UNICOR/unicor/resist_sabah_example_pro.tif"
# Save as tiff
with rio.open(pathout, 'w', **profile) as dst:
    dst.write(data)

# Check that the file was saved as geotiff
with rio.open(pathout, 'r') as src:
    print(src.profile)

# Convert resistance grid to suitability
iR = 'C:/Users/pj276/Projects/UNICOR_arch/unicor/size6_side1000px_1000000totpix.tif'
with rio.open(iR) as src:
    r = src.read(1)
    r[r==-9999] = np.nan
    r = 1/r*100
    r[np.isnan(r)] = -9999
    r = np.expand_dims(r, axis=0)
    profile = src.profile

iRout = 'C:/Users/pj276/Projects/UNICOR_arch/unicor/sabah_test_suitability.tif'
with rio.open(iRout, 'w', **profile) as dst:
    dst.write(r)

with rio.open(iRout) as src:
    r = src.read(1)

# Get all pairwise sums of corridors
# Each sum should then be a minimum conditional transit cost surface
# Thresholding should then get you a corridor
# This works but not sure if it scales
import numpy as np
from itertools import combinations
a = np.array([[1,2,3],[4,5,6],[7,8,9]])
c = ([a[:,i]+a[:,j] for i,j in combinations(range(a.shape[1]),2)])



# CHeck raster array
rname = 'C:/Users/pj276/Projects/UNICOR_arch/unicor/roads.tif'
with rio.open(rname, 'r') as src:
    profile = src.profile
    rio._io.Statistics(src,'min')
    data = src.read()

profile.update(nodata=-9999)

iRout = 'C:/Users/pj276/Projects/UNICOR_arch/unicor/roads_ndmod.tif'
with rio.open(iRout, 'w', **profile) as dst:
    dst.write(data)

# Data frame with coordinates to shapefile
# Convert to pandas df
df = pd.read_csv('C:/Users/pj276/OneDrive - Northern Arizona University/All/Proposals/NASA_GEDI_2023/Amazon_series-AMZ_CAMTRAP_v1.0/LEEClab-Amazon_series-6f7b14d/Amazon_camtrap/DATASET_09_02_2022/AMZ_CAMTRAP_UNIT.csv', encoding = 'ISO-8859-1')
df = df.loc[df.COL_END_YR >= 2017]
df = df.loc[df.BAIT == 0]
df = df[~df['VEG_LANDUSE_TYPE_POINT'].str.contains('Urban area', na=False)]
df = df.loc[(df.COUNTRY == 'Brazil') | (df.COUNTRY == 'France')]
#df = df.loc[df.COUNTRY == 'France']
df = df[~df['DATA_TEAM'].str.contains('Schaub', na=False)]
df = df[~df['ORDEM_BD'].str.contains('ANDCO_', na=False)]
df = df[~df['ORDEM_BD'].str.contains('AGUSIL_', na=False)]
df = df[~df['RECORD_ID'].str.contains('Paleosuchus_', na=False)]


# Convert to geodataframe
gdf = gpd.GeoDataFrame(df, geometry=gpd.points_from_xy(df.LONG_X_POINT, df.LAT_Y_POINT), crs="epsg:4326")
ofile = 'C:/Users/pj276/OneDrive - Northern Arizona University/All/Proposals/NASA_GEDI_2023/Amazon_series-AMZ_CAMTRAP_v1.0/LEEClab-Amazon_series-6f7b14d/Amazon_camtrap/DATASET_09_02_2022/AMZ_CAMTRAP_UNIT_filter_v2.shp'
gdf.to_file(ofile)

# Cumulative sum by year
df = pd.read_csv('C:/Users/pj276/OneDrive - Northern Arizona University/All/Proposals/NASA_GEDI_2023/Amazon_series-AMZ_CAMTRAP_v1.0/LEEClab-Amazon_series-6f7b14d/Amazon_camtrap/DATASET_09_02_2022/AMZ_CAMTRAP_UNIT.csv', encoding = 'ISO-8859-1')
df = df.loc[df.BAIT == 0]
df = df[~df['VEG_LANDUSE_TYPE_POINT'].str.contains('Urban area', na=False)]
df = df[~df['DATA_TEAM'].str.contains('Schaub', na=False)]
#df = df.loc[(df.COUNTRY == 'Brazil') | (df.COUNTRY == 'France')]
df = df.loc[ df.COUNTRY == 'France']
df = df.groupby('COL_END_YR')['COL_END_YR'].count()
df.plot()
print (df)

# Read in csv of lat/lon of paisagen lidar data
df = pd.read_csv(r"C:\Users\pj276\OneDrive - Northern Arizona University\All\Proposals\NASA_GEDI_2023\cms_brazil_lidar_tile_inventory.csv")
df["geometry"] = df.apply(lambda row: 'POLYGON(('+str(row.min_lon)+' '+str(row.max_lat)+','+str(row.max_lon)+' '+str(row.max_lat)+','+str(row.max_lon)+' '+str(row.min_lat)+','+str(row.min_lon)+' '+str(row.min_lat)+','+str(row.min_lon)+' '+str(row.max_lat)+'))', axis=1)
df['geometry'] = gpd.GeoSeries.from_wkt(df['geometry'])
geo_df = gpd.GeoDataFrame(df, geometry='geometry',crs='epsg:4326')
ofile = 'C:/Users/pj276/OneDrive - Northern Arizona University/All/Proposals/NASA_GEDI_2023/paisagens_bboxes.shp'
geo_df.to_file(ofile)


from shapely.geometry import shape
polys = []
for shape, value in shapes(segments, transform=affine):
    polys.append(shape)

geom = [shape(i) for i in polys]
gpd.GeoDataFrame({'geometry':geom})
                         geometry

ss = gpd.read_file(r"C:\Users\pj276\OneDrive - Northern Arizona University\All\Proposals\NASA_GEDI_2023\SnapshotSerengetiBboxes_20190903\SnapshotSerengetiBboxes_20190903.json")
ss = pd.read_json(r"C:\Users\pj276\OneDrive - Northern Arizona University\All\Proposals\NASA_GEDI_2023\SnapshotSerengetiBboxes_20190903\SnapshotSerengetiBboxes_20190903.json")