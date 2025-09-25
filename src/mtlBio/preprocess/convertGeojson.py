
from pathlib import Path
import duckdb
from shapely import wkb, wkt
from shapely.geometry import Point
import geopandas as gpd
from mtlBio.config import configs

def convertParquetGeojson(file_path:str = None, merge_col:str = None, debug:bool = False):
    print("Geographic data, exporting to GeoJSON")

    file_path = Path(file_path)
    output_path = Path(f"{file_path.parent}/{file_path.stem}.geojson")
    shp_file_path = Path(f"{configs.data_dir}/interim/geospatial/{file_path.stem}.shp")
    
    con = duckdb.connect()
    con.execute("INSTALL spatial;")
    con.execute("LOAD spatial;")
    
    df = con.execute(f""" SELECT * FROM '{file_path}'""").df()    
    gdf = gpd.read_file(shp_file_path)
    gdf = gdf.merge(df, on = merge_col)
    
    #Clean Nan values 
    gdf = gdf.dropna(how="any").reset_index(drop=True)

    try:
        gdf.to_file(output_path, driver = 'GeoJSON')
        print(f'Successfully saved {output_path}')
    except Exception as e:
        print(f'Failed to save {output_path}: {e}')

if __name__ == "__main__":
    
    quartiers_file = f"{configs.data_dir}/processed/quartiers.parquet"
    grid_file = f"{configs.data_dir}/processed/grid.parquet"
    convertParquetGeojson(quartiers_file, "qrt_id")
    convertParquetGeojson(grid_file, "grid_id")