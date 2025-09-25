
from pathlib import Path
import duckdb
from shapely import wkb, wkt
from shapely.geometry import Point
import geopandas as gpd
from mtlBio.config import configs

def convertParquetGeojson(file_path:str = None, debug:bool = False, merge_col:str = None):
    print("Geographic data, exporting to GeoJSON")

    file_path = Path(file_path)
    output_path = Path(f"{file_path.parent}/{file_path.stem}.geojson")
    shp_file_path = Path(f"{configs.data_dir}/interim/geospatial/{file_path.stem}.shp")
    
    con = duckdb.connect()
    con.execute("INSTALL spatial;")
    con.execute("LOAD spatial;")
    
    df = con.execute(f""" SELECT * FROM '{file_path}'""").df()

    print(df)
    
    gdf = gpd.read_file(shp_file_path)

    print(gdf)
    
    gdf = gdf.merge(df, on = 'qrt_id')
    print(gdf)
    
    #Clean Nan values 
    gdf = gdf.dropna(how="any").reset_index(drop=True)

    try:
        gdf.to_file(output_path, driver = 'GeoJSON')
        print(f'Successfully saved {output_path}')
    except Exception as e:
        print(f'Failed to save {output_path}: {e}')

if __name__ == "__main__":
    
    file = f"{configs.data_dir}/processed/quartiers.parquet"
    convertParquetGeojson(file)