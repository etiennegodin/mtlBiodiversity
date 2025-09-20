import streamlit_folium
import streamlit as st
import folium
import geopandas as gpd
from pathlib import Path
from folium.features import GeoJsonTooltip, GeoJsonPopup


gjson_path = Path("data/processed/parks.geojson")

gdf = gpd.read_file(gjson_path)

m = folium.Map(location=[45.5017, -73.5673], zoom_start=11) 
choropleth = folium.Choropleth(
    geo_data=gdf.__geo_interface__,
    data=gdf,
    columns=["park_name","species_richness",],
    key_on="feature.properties.park_name",
    fill_color="viridis",
    fill_opacity=0.7,
    line_opacity=0.2,
    legend_name="Species Richness"
).add_to(m)

# Add hover tooltips
folium.GeoJsonTooltip(
    fields=["park_name", "species_richness"],
    aliases=["Park:", "Species Richness:"],
).add_to(choropleth.geojson)


streamlit_folium.st_folium(m, width=700, height=500)