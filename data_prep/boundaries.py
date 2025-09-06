import geopandas as gpd

park_file = "data/raw/espace_verts/Espace_Vert.shp"

gdf = gpd.read_file(park_file)
