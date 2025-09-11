import pandas as pd
import geopandas as gpd
from pathlib import Path
import streamlit as st 


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



gdf = gpd.read_file(gjson_path)

# Side bar
with st.sidebar:
    st.title = TITLE
    park_list = list(gdf.park_name.unique().tolist())

    selected_park = st.selectbox('Select park', park_list, index = len(park_list)-1)
    color_theme_list = ['blues', 'cividis', 'greens', 'inferno', 'magma', 'plasma', 'reds', 'rainbow', 'turbo', 'viridis']

    selected_color_theme = st.selectbox('Select a color theme', color_theme_list)

# Layout 
col = st.columns((1.4,4.5,2), gap = 'medium')


