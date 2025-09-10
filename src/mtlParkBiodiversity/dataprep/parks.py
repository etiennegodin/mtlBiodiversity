import geopandas as gpd
from pathlib import Path
from ..core import clip_to_region, convert_crs
from . import target_crs

raw_data_path = Path("data/raw" )
output_path = Path("data/interim/parks")

Path.mkdir(output_path, exist_ok= True)


#Read boundary file 


def save_gdf(gdf, output_file_path):
    try:
        gdf.to_file(output_file_path)
        print(f'Successfuly saved {output_file_path}')
        return True
    except Exception as e:
        print(f"Failed to save {output_file_path.stem} : {e}")
        return False

def process_shp(file :Path = None, output_file_path :Path = None, city_boundary_gdf :gpd.GeoDataFrame = None):

        # Read shapefile 
    gdf = gpd.read_file(file)

    #Clip gdf 
    gdf_clipped = clip_to_region(gdf, city_boundary_gdf)

    # Convert to target crs before saving 
    gdf_out = convert_crs(gdf_clipped, target_crs = target_crs)
    if not save_gdf(gdf_out,output_file_path):
        return False
    return True

def process_gpgk(file :Path = None, output_file_path :Path = None, city_boundary_gdf :gpd.GeoDataFrame = None, layer = None):
    #Create new filename base on layer name
    output_file_name = file.stem + f'_clipped_{layer}.shp'
    # Final output path with new subfolder
    output_file_path = output_file_path / output_file_name

    #Read gpgk layer 
    gdf = gpd.read_file(file, layer = layer)
    #Clip gdf 
    gdf_clipped = clip_to_region(gdf, city_boundary_gdf)

    # Convert to target crs before saving 
    gdf_out = convert_crs(gdf_clipped, target_crs = target_crs)
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


def prep_parks(force = False):

    city_boundary_file = [f for f in raw_data_path.rglob("*.shp")][0]  # Assuming there's only one .shp file for the city boundary
    city_boundary_gdf = gpd.read_file(city_boundary_file)

    parks_path = raw_data_path / 'parks'

    parks_files = [f for f in parks_path.rglob("*") if f.suffix in ['.shp', '.gpkg']]

    print(parks_files)


    for file in parks_files:
        # Set target output_path 
        output_file_name = file.stem + '_clipped' + file.suffix 
        output_file_path = output_path / output_file_name

        if output_file_path.exists and force:

            if file.suffix == '.gpkg':

                new_output_path = output_path / file.stem  # If gpgk, modify output in a subfolder for better clarity 
                new_output_path.mkdir(exist_ok= True) # Create subfolder if not existing
                
                #Create list of gpgk layers
                layers = list_gpgk_layers(file)
                #Iterate over each layer and export clipped shp 
                for idx, layer in enumerate(layers):
                    print(layer, f" : {idx+1} / {len(layers)}")

                    process_gpgk(file, layer, output_file_path = new_output_path)

            elif file.suffix == '.shp':

                process_shp(file, output_file_path, city_boundary_gdf = city_boundary_gdf)

        else:
            print(f"# File {output_file_path} already exists, skipping.")
    
    return True

if __name__ == "__main__":
    prep_parks(force= False)