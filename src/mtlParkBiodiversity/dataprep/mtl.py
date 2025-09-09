import geopandas as gpd
from pathlib import Path
from ..core import clip_to_region, convert_crs
from ..dataprep import target_crs

raw_data_path = Path("data/raw" )
output_path = Path("data/interim")

foretOuverte_file = raw_data_path / "PRODUITS_IEQM_31H_GPKG"/ "PRODUITS_IEQM_31H.gpkg"

mtl_boundary_file = raw_data_path / "limites-administratives-agglomeration-nad83" / "limites-administratives-agglomeration-nad83.shp"
park_file_path = raw_data_path / "espace_verts" / "Espace_Vert.shp"
bois_file_path = raw_data_path / "bois"/ "bois.shp"

mtl_data = [park_file_path, bois_file_path]

#Read boundary file 

mtl_boundary_gdf = gpd.read_file(mtl_boundary_file)

def save_gdf(gdf, output_file_path):
    try:
        gdf.to_file(output_file_path)
        print(f'Successfuly saved {output_file_path}')
        return True
    except Exception as e:
        print(f"Failed to save {output_file_path.stem} : {e}")
        return False
    
def extract_gpgk_layers(file : Path, slice : tuple = None):

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


def prep_mtl(force = False):

    for file in mtl_data:
        # Set target output_path 
        output_file_name = file.stem + '_clipped' + file.suffix 
        output_file_path = output_path / output_file_name

        if output_file_path.exists and force:

            if file.suffix == '.gpkg':

                new_output_path = output_path / file.stem  # If gpgk, modify output in a subfolder for better clarity 
                new_output_path.mkdir(exist_ok= True) # Create subfolder if not existing
                
                #Create list of gpgk layers
                layers = extract_gpgk_layers(file)
                #Iterate over each layer and export clipped shp 
                for idx, layer in enumerate(layers):
                    print(layer, f" : {idx+1} / {len(layers)}")

                    #Create new filename base on layer name
                    output_file_name = file.stem + f'_clipped_{layer}.shp'
                    # Final output path with new subfolder
                    output_file_path = new_output_path / output_file_name

                    #Read gpgk layer 
                    gdf = gpd.read_file(file, layer = layer)
                    #Clip gdf 
                    gdf_clipped = clip_to_region(gdf, mtl_boundary_gdf)

                    # Convert to target crs before saving 
                    gdf_out = convert_crs(gdf_clipped, target_crs = target_crs)

                    if not save_gdf(gdf_out,output_file_path):
                        continue

            elif file.suffix == '.shp':

                # Read shapefile 
                gdf = gpd.read_file(file)

                #Clip gdf 
                gdf_clipped = clip_to_region(gdf, mtl_boundary_gdf)

                # Convert to target crs before saving 
                gdf_out = convert_crs(gdf_clipped, target_crs = target_crs)
                if not save_gdf(gdf_out,output_file_path):
                    continue

        else:
            print(f"# File {output_file_path} already exists, skipping.")

if __name__ == "__main__":
    prep_mtl(force= False)