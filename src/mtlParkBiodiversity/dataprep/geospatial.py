import geopandas as gpd
from pathlib import Path
from ..core import convert_crs
from ..column_mapper import unify_columns
from . import target_crs

RAW_DATA_PATH = Path("data/raw" )
OUTPUT_PATH = Path("data/interim/geospatial")

GEOSPATIAL_DATA_COL = ['geospatial_id', 'geospatial_name', 'geospatial_type1', 'geometry']

#Read boundary file 
def save_gdf(gdf, output_file_path):
    try:
        gdf.to_file(output_file_path)
        print(f'Successfuly saved {output_file_path}')
        return True
    except Exception as e:
        print(f"Failed to save {output_file_path.stem} : {e}")
        return False

def process_shp(file :Path = None, output_file_path :Path = None, force = False):

    # Read shapefile 
    gdf = gpd.read_file(file)
    
    # Format columns 
    gdf, columns = unify_columns(file, force = force)

    #Add geometry back if missing
    if "geometry" not in columns:
        columns.append('geometry')

    # Keep only mapped columns
    gdf = gdf[columns]

    print('Geospatial data processed, saving')
    # Convert to target crs before saving 
    gdf_out = convert_crs(gdf, target_crs = target_crs)
    if not save_gdf(gdf_out,output_file_path):
        return False
    return True

def process_gpgk(file :Path = None, output_file_path :Path = None, layer = None):
    #Create new filename base on layer name
    output_file_name = file.stem + f'_clipped_{layer}.shp'
    # Final output path with new subfolder
    output_file_path = output_file_path / output_file_name

    #Read gpgk layer 
    gdf = gpd.read_file(file, layer = layer)

    # Convert to target crs before saving 
    gdf_out = convert_crs(gdf, target_crs = target_crs)
    if not save_gdf(gdf_out,output_file_path):
        return False
    return True

def list_gpgk_layers(file : Path, slice : tuple = None):

    # Create df from layers of gpgk 
    layers_df = gpd.list_layers(file)

    #Keep only geometry layers 
    values_to_remove = ['MultiPolygon', 'Polygon']
    df_filtered = layers_df[layers_df['geometry_type'].isin(values_to_remove)]
    layers = df_filtered['name'].tolist()

    #Splice list if specified
    if slice:
        layers = layers[slice[0]:slice[1]]

    return layers


def prep_geospatial(force = False):

    # Create output directory if not existing
    Path.mkdir(OUTPUT_PATH, exist_ok= True)

    # List all park files (shp, gpgk) in RAW_DATA_PATH / spatial
    geospatial_path = RAW_DATA_PATH / 'geospatial'
    geospatial_files = [f for f in geospatial_path.rglob("*") if f.suffix in ['.shp', '.gpkg']]

    # Iterate over each park file and process
    for file in geospatial_files:
        # Set target OUTPUT_PATH
        output_file_name = f'{file.stem}_inter.shp' 
        output_file_path = OUTPUT_PATH / output_file_name

        if (output_file_path.exists() and force) or (not output_file_path.exists()):

            if file.suffix == '.shp':
                process_shp(file, output_file_path, force = force)

            """
            elif file.suffix == '.gpkg':

                new_OUTPUT_PATH = OUTPUT_PATH / file.stem  # If gpgk, modify output in a subfolder for better clarity 
                new_OUTPUT_PATH.mkdir(exist_ok= True) # Create subfolder if not existing
                
                #Create list of gpgk layers
                layers = list_gpgk_layers(file)
                #Iterate over each layer and export clipped shp 
                for idx, layer in enumerate(layers):
                    print(layer, f" : {idx+1} / {len(layers)}")

                    process_gpgk(file, layer, output_file_path = new_OUTPUT_PATH)
            """
        else:
            print(f"# File {output_file_path} already exists, skipping.")
    
    return True

if __name__ == "__main__":
    prep_geospatial(force= False)