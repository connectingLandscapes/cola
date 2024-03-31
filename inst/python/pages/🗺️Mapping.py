#streamlit version of the geemap, developed with folium
import streamlit as st
import ee
import folium
from streamlit_folium import folium_static

# Initialize the Earth Engine module.
ee.Initialize()

# Myanmar kernel and lcp layers. Mask out zeros.
# Todo -- change the source of kernel and lcp layers
baseCore = ee.Image('')
baseCore = baseCore.updateMask(baseCore.gt(0))
baseLCP = ee.Image('')
baseLCP = baseLCP.updateMask(baseLCP.gt(0))
growthCore = ee.Image('')
growthCore = growthCore.updateMask(growthCore.gt(0))
growthLCP = ee.Image('')
growthLCP = growthLCP.updateMask(growthLCP.gt(0))
baseSuitability = ee.Image('')
baseSuitability = baseSuitability.updateMask(baseSuitability.gt(0))

# Ecoregions clipped to Colombia boundary
aoi = ee.FeatureCollection("projects/ee-jantzenator/assets/gedi_col_ecoregion_sset/colombia_clip") 

# Country boundary
aoi = ee.FeatureCollection("FAO/GAUL/2015/level0")
# Filter
mmr = aoi.filter(ee.Filter.eq('ADM0_NAME', 'Myanmar'))

# WDPA
pa = ee.FeatureCollection("WCMC/WDPA/current/polygons")
pa = pa.filter(ee.Filter.eq('ISO3', 'MMR'))
# Could add IUCN_CAT as well (I, II, III, etc.)
# Create an empty image into which to paint the features, cast to byte.
empty = ee.Image().byte()
pap = empty.paint({
  'featureCollection': pa,
  'color': 1
})


# Create a folium map object.
my_map = folium.Map(location=[21.9162, 95.9560], zoom_start=5)  # Centered on Myanmar

# Define color palettes as lists of hex color strings.
viridis_core = ['#481567FF', '#482677FF', '#453781FF', '#404788FF', '#39568CFF',
                '#33638DFF', '#2D708EFF', '#287D8EFF', '#238A8DFF', '#1F968BFF',
                '#20A387FF', '#29AF7FFF', '#3CBB75FF', '#55C667FF', '#73D055FF',
                '#95D840FF', '#B8DE29FF', '#DCE319FF', '#FDE725FF']

# For now, let's assume all color palettes are the same
viridis_lcp = viridis_suitability = viridis_core

# Add checkboxes to the sidebar for layer visibility.
pap_check = st.sidebar.checkbox('Protected Areas', value=False)
base_suitability_check = st.sidebar.checkbox('Habitat Suitability', value=False)
base_core_check = st.sidebar.checkbox('Core Habitat (baseline)', value=False)
growth_core_check = st.sidebar.checkbox('Core Habitat (growth)', value=False)
base_lcp_check = st.sidebar.checkbox('Connectivity (baseline)', value=False)
growth_lcp_check = st.sidebar.checkbox('Connectivity (growth)', value=False)

# Add sliders to the sidebar for layer opacity.
pap_opacity = st.sidebar.slider('Opacity for Protected Areas', min_value=0.0, max_value=1.0, value=1.0, step=0.01)
base_suitability_opacity = st.sidebar.slider('Opacity for Habitat Suitability', min_value=0.0, max_value=1.0, value=1.0, step=0.01)
base_core_opacity = st.sidebar.slider('Opacity for Core Habitat (baseline)', min_value=0.0, max_value=1.0, value=1.0, step=0.01)
growth_core_opacity = st.sidebar.slider('Opacity for Core Habitat (growth)', min_value=0.0, max_value=1.0, value=1.0, step=0.01)
base_lcp_opacity = st.sidebar.slider('Opacity for Connectivity (baseline)', min_value=0.0, max_value=1.0, value=1.0, step=0.01)
growth_lcp_opacity = st.sidebar.slider('Opacity for Connectivity (growth)', min_value=0.0, max_value=1.0, value=1.0, step=0.01)

# If the checkbox is checked, add the corresponding layer to the map.

if pap_check:
    # Add the 'pap' layer to the map with the specified opacity.
    # This is a placeholder. You'll need to replace it with the actual code to add the layer.
    my_map.add_layer(pap, opacity=pap_opacity)

if base_suitability_check:
    # Add the 'baseSuitability' layer to the map with the specified opacity.
    # This is a placeholder. You'll need to replace it with the actual code to add the layer.
    my_map.add_layer(baseSuitability, opacity=base_suitability_opacity)

if base_core_check:
    # Add the 'baseCore' layer to the map with the specified opacity.
    # This is a placeholder. You'll need to replace it with the actual code to add the layer.
    my_map.add_layer(baseCore, opacity=base_core_opacity)

if growth_core_check:
    # Add the 'growthCore' layer to the map with the specified opacity.
    # This is a placeholder. You'll need to replace it with the actual code to add the layer.
    my_map.add_layer(growthCore, opacity=growth_core_opacity)

if base_lcp_check:
    # Add the 'baseLCP' layer to the map with the specified opacity.
    # This is a placeholder. You'll need to replace it with the actual code to add the layer.
    my_map.add_layer(baseLCP, opacity=base_lcp_opacity)

if growth_lcp_check:
    # Add the 'growthLCP' layer to the map with the specified opacity.
    # This is a placeholder. You'll need to replace it with the actual code to add the layer.
    my_map.add_layer(growthLCP, opacity=growth_lcp_opacity)

# Create a folium map
m = folium.Map(location=[45.5236, -122.6750])

# Display the map in Streamlit
folium_static(m)