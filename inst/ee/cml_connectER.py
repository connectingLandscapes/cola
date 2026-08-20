"""
EarthRanger → CoLa DSS  |  Live Data Integration Demo
======================================================
This script demonstrates how the Connecting Landscapes (CoLa) DSS could pull
animal tracking data directly from EarthRanger (via the Ecoscope Python library)
instead of requiring users to manually upload files.

The two outputs — a habitat suitability raster (GeoTIFF) and a source points
shapefile — are exactly the inputs CoLa's DSS expects for connectivity modelling.

Dependencies
------------
    pip install ecoscope geopandas rasterio scipy numpy

EarthRanger is a wildlife monitoring platform used by 100+ conservation
organisations globally. Ecoscope is its official open-source Python API:
    https://github.com/wildlife-dynamics/ecoscope

Author:  Giraffe Conservation Foundation  |  KAZA TFCA DSS Connectivity Project
"""

import os
import numpy as np
import geopandas as gpd
import rasterio
from rasterio.transform import from_origin
from scipy.stats import gaussian_kde
from shapely.geometry import Point
from ecoscope.io.earthranger import EarthRangerIO
import ime, os, sys, argparse
from datetime import datetime

# =============================================================================
# 1. CONNECTION SETTINGS
#    EarthRanger is a cloud platform — each conservation organisation has its
#    own server. No data upload needed; we query it live.
# =============================================================================
#
# Usage:
#    python crk_scipy.py <server> <username> <pwd> <datefrom> <dateto>
#                        <subject> <outdir> <res> <maxLive> <upCRS>
# Arguments:
#    server        : Organisation server URL
#    username      : Username
#    pwd           : Password
#    datefrom      : Date-from, format YYYY-mm-dd
#    dateto        : Date-to, format YYYY-mm-dd
#    subject       : The subject group to analyse (species + study area, as defined in EarthRanger)
#    outdir        : Output dir
#    res           : Raster resolution


def parse_args():
    p = argparse.ArgumentParser(description='SDM MODIS Wall-to-Wall Prediction')
    # ── Core identifiers (shared across scripts) ──────────────────────────────
    p.add_argument('--server', type = str, required = True,
        help='Folder in your Google Drive account to save the results'
    p.add_argument('--username',      type=str, required=True)
    p.add_argument('--pwd',      type=str, required=True,
                   help='Your password')
    p.add_argument('--datefrom', type=str, required=True,
                   help='Initital date, format YYYY-mm-dd')
   p.add_argument('--dateto', type=str,  required=True,
                   help='Final date, format YYYY-mm-dd')
    p.add_argument('--subjetgroup',  type=str, required=True,
                    help = 'The subject group to analyse (species + study area, as defined in EarthRanger)')
    p.add_argument('--oudir',        type=str required=True,
                    help = 'Local path where results will be saved')
    p.add_argument('--spatresinmeters',   type=int, required=True,
                   help = 'KDE raster resolution in metres')
    p.add_argument('--region',  type=str, required=True,
                   help='Feature with the extent to export the layers')
    return p.parse_args()

args = parse_args()
ER_SERVER = args.server
ER_USERNAME = args.username
ER_PASSWORD = args.pwd
DATE_FROM = args.datefrom
DATE_TO = args.dateto
SUBJECT_GROUP = args.subjetgroup
OUTPUT_DIR = args.oudir
KDE_RES_M = args.spatresinmeters

# ER_SERVER   = sys.argv[1] # "https://your-organisation.pamdas.org"   # ← organisation server URL
# ER_USERNAME = sys.argv[2] # "your_username"
# ER_PASSWORD = sys.argv[3] # "your_password"
## Date range for the analysis
#DATE_FROM = sys.argv[4] # "2025-01-01"
#DATE_TO   = sys.argv[5] # "2025-04-01"
## The subject group to analyse (species + study area, as defined in EarthRanger)
#SUBJECT_GROUP = sys.argv[6] # "Your_Species_StudyArea"
## Output folder for CoLa inputs
#OUTPUT_DIR    = sys.argv[7] # "./cola_inputs"
#KDE_RES_M     = sys.argv[8] # 500       # KDE raster resolution in metres

# =============================================================================
# 2. CONNECT AND PULL DATA
#    EarthRangerIO handles authentication, pagination, and rate limiting.
#    A single call returns all relocations for the chosen group and date range.
# =============================================================================

def pull_tracking_data(server, username, password, group, date_from, date_to):
    """Connect to EarthRanger and return a GeoDataFrame of GPS relocations."""

    print(f"Connecting to EarthRanger: {server}")
    er = EarthRangerIO(server=server, username=username, password=password)

    print(f"Pulling relocations  |  group: {group}  |  {date_from} → {date_to}")
    gdf = er.get_subjectgroup_observations(
        subject_group_name = group,
        since              = f"{date_from}T00:00:00Z",
        until              = f"{date_to}T23:59:59Z",
        relocations        = False,        # return plain GeoDataFrame
    )

    if gdf is None or len(gdf) == 0:
        raise ValueError("No relocations returned. Check group name and date range.")

    print(f"  Retrieved {len(gdf):,} relocations.")
    return _ensure_wgs84(gdf)


def _ensure_wgs84(gdf):
    """Make sure the GeoDataFrame has valid WGS84 point geometry."""
    if "geometry" not in gdf.columns or gdf.geometry.is_empty.all():
        lon = next((c for c in ["longitude", "lon", "lng"] if c in gdf.columns), None)
        lat = next((c for c in ["latitude",  "lat"]        if c in gdf.columns), None)
        if not (lon and lat):
            raise ValueError(f"No geometry or coordinate columns found. Available: {list(gdf.columns)}")
        gdf["geometry"] = [Point(r[lon], r[lat]) for _, r in gdf.iterrows()]
        gdf = gpd.GeoDataFrame(gdf, geometry="geometry", crs="EPSG:4326")
    if gdf.crs is None:
        gdf = gdf.set_crs("EPSG:4326")
    elif gdf.crs.to_epsg() != 4326:
        gdf = gdf.to_crs("EPSG:4326")
    return gdf


# =============================================================================
# 3. REPROJECT TO UTM
#    CoLa requires a metric projected CRS (not lat/lon).
#    We auto-detect the appropriate UTM zone from the data centroid.
# =============================================================================

def reproject_utm(gdf):
    """Reproject a WGS84 GeoDataFrame to the appropriate UTM zone."""
    centroid = gdf.geometry.unary_union.centroid
    lon, lat  = centroid.x, centroid.y
    zone      = int((lon + 180) / 6) + 1
    epsg      = f"EPSG:{'326' if lat >= 0 else '327'}{zone:02d}"
    print(f"  Reprojecting to {epsg}")
    return gdf.to_crs(epsg), epsg


# =============================================================================
# 4. EXPORT: SOURCE POINTS SHAPEFILE
#    CoLa uses these as the starting locations for least-cost path analysis.
# =============================================================================

def export_points(gdf, out_dir, name):
    """Save GPS locations as a shapefile with X/Y attribute columns."""
    path = os.path.join(out_dir, f"{name}_points.shp")

    pts       = gdf.copy()
    pts["X"]  = pts.geometry.x
    pts["Y"]  = pts.geometry.y
    keep      = ["X", "Y", "geometry"]
    for col in ["subject_name", "recorded_at"]:
        if col in pts.columns:
            keep.append(col)

    export         = pts[keep].copy()
    export.columns = [c[:10] for c in export.columns]   # shapefile 10-char limit
    export.to_file(path)

    print(f"  Source points → {path}  ({len(export):,} points)")
    return path


# =============================================================================
# 5. EXPORT: HABITAT SUITABILITY RASTER (KDE)
#    Kernel Density Estimation converts GPS point density into a continuous
#    0–1 habitat suitability surface — the standard input for CoLa's
#    resistance modelling step.
#
#    IMPORTANT: CoLa requires square pixels. We use from_origin() with a fixed
#    resolution to guarantee identical X and Y cell dimensions.
# =============================================================================

def export_kde_raster(gdf, epsg, out_dir, name, resolution_m=500):
    """Compute a KDE density surface and save as a GeoTIFF."""
    path = os.path.join(out_dir, f"{name}_hs_kde.tif")

    xs = gdf.geometry.x.values.astype(float)
    ys = gdf.geometry.y.values.astype(float)
    mask = np.isfinite(xs) & np.isfinite(ys)
    xs, ys = xs[mask], ys[mask]

    if len(xs) < 5:
        raise ValueError("Too few valid coordinates for KDE (minimum 5 required).")

    # Buffered extent so animals near the edge are fully captured
    buf    = max((xs.max() - xs.min()), (ys.max() - ys.min())) * 0.10 + resolution_m
    xmin, xmax = xs.min() - buf, xs.max() + buf
    ymin, ymax = ys.min() - buf, ys.max() + buf

    # Grid dimensions — fixed resolution guarantees square pixels
    ncols = max(100, min(2000, int(np.ceil((xmax - xmin) / resolution_m))))
    nrows = max(100, min(2000, int(np.ceil((ymax - ymin) / resolution_m))))

    grid_x  = xmin + (np.arange(ncols) + 0.5) * resolution_m
    grid_y  = ymin + (np.arange(nrows) + 0.5) * resolution_m
    xx, yy  = np.meshgrid(grid_x, grid_y)

    kde     = gaussian_kde(np.vstack([xs, ys]))
    density = kde(np.vstack([xx.ravel(), yy.ravel()])).reshape(nrows, ncols).astype("float32")

    # Normalise to 0–1 (CoLa habitat suitability range)
    d_min, d_max = density.min(), density.max()
    if d_max > d_min:
        density = (density - d_min) / (d_max - d_min)

    density = np.flipud(density)   # rasterio origin is top-left

    with rasterio.open(
        path, "w",
        driver    = "GTiff",
        height    = nrows,
        width     = ncols,
        count     = 1,
        dtype     = "float32",
        crs       = epsg,
        transform = from_origin(xmin, ymin + nrows * resolution_m, resolution_m, resolution_m),
        nodata    = -9999.0,
    ) as dst:
        dst.write(density, 1)

    print(f"  Habitat suitability raster → {path}  ({ncols}×{nrows} px at {resolution_m} m)")
    return path


# =============================================================================
# 6. RUN
# =============================================================================

def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # Sanitise group name for use as a filename
    safe_name = "".join(c if c.isalnum() or c in "_-" else "_" for c in SUBJECT_GROUP)

    print("\n" + "═" * 60)
    print("  EarthRanger → CoLa  |  Live Data Integration Demo")
    print("═" * 60 + "\n")

    # Pull live data from EarthRanger
    gdf          = pull_tracking_data(ER_SERVER, ER_USERNAME, ER_PASSWORD,
                                      SUBJECT_GROUP, DATE_FROM, DATE_TO)
    gdf_utm, epsg = reproject_utm(gdf)

    print()

    # Export CoLa-ready inputs
    shp_path = export_points(gdf_utm,              OUTPUT_DIR, safe_name)
    tif_path = export_kde_raster(gdf_utm, epsg,    OUTPUT_DIR, safe_name, KDE_RES_M)

    print(f"""
{"═" * 60}
  Done. CoLa inputs ready:

  Habitat suitability raster : {tif_path}
  Source points shapefile    : {shp_path}

  These files can be loaded directly into the CoLa DSS
  (Habitat suitability ↔ resistance  →  Kernels / Corridors).

  Integration note for developers:
  ─────────────────────────────────
  To embed this in the DSS, the user would supply:
    • EarthRanger server URL  (their organisation's instance)
    • Username / password     (or OAuth token)
    • Subject group name      (dropdown from er.get_subjectgroups())
    • Date range
    • KDE resolution

  The DSS would call pull_tracking_data() and both export
  functions, then pass the outputs directly to the existing
  resistance / kernel / corridor pipeline — no file upload needed.
{"═" * 60}
""")


if __name__ == "__main__":
    main()
