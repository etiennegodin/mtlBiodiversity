import geopandas as gpd
from pathlib import Path
import json
from mtlBio.core import convert_crs, convertToPath

#Read boundary file 
def save_gdf(gdf, output_file_path):
    try:
        gdf.to_file(output_file_path)
        print(f'Successfuly saved {output_file_path}')
        return True
    except Exception as e:
        print(f"Failed to save {output_file_path.stem} : {e}")
        return False

def process_shp(input_file :str = None, output_file :str = None, crs = None):

    input_path = convertToPath(input_file)
    output_path = convertToPath(output_file)

    # Read shapefile 
    gdf = gpd.read_file(input_path)
    
    #Load json 
    json_mapping = f"{input_path.parent}/{input_path.stem}_column_mapping.json"
    with open(json_mapping, 'r') as f:
        mapping = json.load(f)
    print(mapping)
    #Remap col name
    gdf = gdf.rename(mapping)

    print('Geospatial data processed, saving')
    # Convert to target crs before saving 
    gdf_out = convert_crs(gdf, target_crs = crs)
    if not save_gdf(gdf_out,output_path):
        return False
    return True
