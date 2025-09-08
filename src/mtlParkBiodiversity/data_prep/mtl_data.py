import geopandas as gpd
from pathlib import Path
from data_prep_functions import clip_to_region

raw_data_path = Path("data/raw" )
output_path = Path("data/interim")

mtl_boundary_file = raw_data_path / "limites-administratives-agglomeration-nad83" / "limites-administratives-agglomeration-nad83.shp"
park_file_path = raw_data_path / "espace_verts" / "Espace_Vert.shp"
bois_file_path = raw_data_path / "bois"/ "bois.shp"

mtl_data = [park_file_path, bois_file_path ]

# Clip 

#Read boundary file 

mtl_boundary_gdf = gpd.read_file(mtl_boundary_file)

for file in mtl_data:

    output_file_name = file.stem + '_clipped' + file.suffix  
    output_file_path = output_path / output_file_name

    gdf = gpd.read_file(file)

    gdf_clipped = clip_to_region(gdf, mtl_boundary_gdf)
    try:
        gdf_clipped.to_file(output_file_path)
        print(f'Successfuly saved {output_file_path}')
    except Exception as e:
        print(f"Failed to save {file.stem} : {e}")