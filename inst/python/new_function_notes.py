# -*- coding: utf-8 -*-
"""
Created on Fri Aug 11 14:17:54 2023

@author: pj276
"""

# For prioritization step, take inverse of factorial least
# cost paths or least cost corridors. 
# Use CRK patches as markers. 
# Use the watershed algorithm in skimage. Flood from the markers.
# This should result in a segmentation of the least cost
# corridors that "fills" up the corridors from their center
# lines to their edges. Each corridor pixel gets filled
# with the label of the nearest CRK patch. This should result
# in relatively even boundaries between CRK patches but
# will need to see how the boundary pixels behave.
# One other question is how the watershed algorithm
# might behave in the presence of pits in the surface.

# Can use morphological operators to identify CRK patches
# that are connected via corridors and then calculate
# importance of corridors using standard approach.
# Instead of the inverse of the factorial least cost surface,
# cost distance from each CRK patch could be calculated and 
# then masked with the factorial least cost surface. When 
# flooded, I expect this might result in a more isotropic
# allocation of pixels to CRK patches. 

aa = pd.read_csv(r"C:\Users\pj276\OneDrive - Northern Arizona University\All\Proposals\NASA_GEDI_2023\sasia_camera_grid_locations.csv")
gdf = gpd.GeoDataFrame(aa, geometry=gpd.points_from_xy(aa.Longitude, aa.Latitude), crs="EPSG:4326")
gdf.to_file(r"C:\Users\pj276\OneDrive - Northern Arizona University\All\Proposals\NASA_GEDI_2023\sasia_camera_grid_locations.shp")

   