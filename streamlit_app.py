import pandas as pd
import geopandas as gpd
from pathlib import Path

import streamlit as st 
import streamlit_folium
import folium


PARENT_FOLDER = Path(__file__).parent
DATA_FOLDER = PARENT_FOLDER / "data/processed/" 

gjson_path = DATA_FOLDER / "parks.geojson"

TITLE = 'Montreal Biodiversity Dashboard'
# Page configs

st.set_page_config(
page_title= TITLE,
layout = 'wide',
initial_sidebar_state= 'expanded'
)

@st.cache_data(ttl='1d')
def get_gdf_data():

    gdf = gpd.read_file(gjson_path)
    return gdf

gdf = get_gdf_data()


# PAGE

f'''
# {TITLE}
Explore Montreal's hidden biodiversity 
[Link](https://inaturalist.ca/home)

'''

'''
'''

# Side bar
with st.sidebar:
    st.title = TITLE
    park_list = list(gdf.park_name.unique().tolist())

    #selected_park = st.selectbox('Select park', park_list, index = len(park_list)-1)
    color_theme_list = ['blues', 'cividis', 'greens', 'inferno', 'magma', 'plasma', 'reds', 'rainbow', 'turbo', 'viridis']

    selected_color_theme = st.selectbox('Select a color theme', color_theme_list)

    
# Layout 
col = st.columns((1.4,4.5,2), gap = 'medium')

with col[1]:
    selected_park = st.selectbox('Which park do you want to look', park_list)

    st.bar_chart(

        gdf,
        x = 'park_name',
        y = 'species_richness'


    )

    m = folium.Map(location=[45.5017, -73.5673], zoom_start=11)

    choropleth_species_richness = folium.Choropleth(
        geo_data=gdf,
        data=gdf,
        name ="Species Richness",
        columns=["park_name","species_richness",],
        key_on="feature.properties.park_name",
        fill_color="YlGn",
        fill_opacity=0.7,
        line_opacity=0.2,
        legend_name="Species Richness",
        show= False
    ).add_to(m)

    choropleth_shannon_index = folium.Choropleth(
    geo_data=gdf,
    data=gdf,
    name ="Shannon Index",
    columns=["park_name","shannon_index",],
    key_on="feature.properties.park_name",
    fill_color="YlGn",
    fill_opacity=0.7,
    line_opacity=0.2,
    legend_name="Shannon Index",
    show = False

    ).add_to(m)


    # Add hover tooltips
    folium.GeoJsonTooltip(
        fields=["park_name", "species_richness"],
        aliases=["Park:", "Species Richness:"],
    ).add_to(choropleth_species_richness.geojson)

    folium.LayerControl().add_to(m)

    streamlit_folium.st_folium(m, width=700, height=500)

