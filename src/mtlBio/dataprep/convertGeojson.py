
from pathlib import Path
import duckdb
from shapely import wkb, wkt
from shapely.geometry import Point
import geopandas as gpd

def list_to_bytes(lst):
    return bytes(lst)  # no int(x,16), just bytes()

def duckdb_geom_to_shapely(x):
    if x is None:
        return None
    # Convert list of ints to bytes, then parse with shapely
    try:
        return wkb.loads(bytes(x))
    except Exception as e:
        # Sometimes DuckDB adds a prefix; skip the first byte
        return wkb.loads(bytes(x[1:]))


def convertParquetGeojson(file_path:str = None, debug:bool = False):
    print("Geographic data, exporting to GeoJSON")

    file_path = Path(file_path)
    output_path = Path(f"{file_path.parent}/{file_path.stem}.geojson")
    
    con = duckdb.connect()
    con.execute("INSTALL spatial;")
    con.execute("LOAD spatial;")
    
    df = con.execute(f""" SELECT * FROM '{file_path}'""").df()
    
    df['geometry'] = df['geom'].apply(lambda lst: wkb.loads(list_to_bytes(lst)))


    #Convert park_geom from wkb to shapely geometry
    #df["geometry"] = df["geom"].apply(lambda x: wkb.loads(bytes(x)))
    #df['geometry'] = df['geom_wkb'].apply(lambda x: wkb.loads(x) if x is not None else None)
    #df['geometry'] = df['geom'].apply(duckdb_geom_to_shapely)
    try:
        df['geometry'] = df['geom'].apply(wkt.loads)
    except:
        df['geometry'] = df['geom'].apply(lambda x: wkb.loads(bytes(x)) if x is not None else None)

        #df['geometry'] = df['geom'].apply(lambda xy: Point(xy))

    print(df)



    #df['geometry'] = df['geom'].apply(lambda x: wkb.loads(bytes(x)) if x is not None else None)

    gdf = gpd.GeoDataFrame(df, geometry = 'geometry', crs = 'EPSG:4326')
    gdf = gdf.drop(columns = ['geom'])
    
    #Clean Nan values 
    gdf = gdf.dropna(how="any").reset_index(drop=True)

    try:
        gdf.to_file(output_path, driver = 'GeoJSON')
        print(f'Successfully saved {output_path}')
    except Exception as e:
        print(f'Failed to save {output_path}: {e}')

if __name__ == "__main__":
    file = "C:/Users/manat/Documents/Projects/mtlBiodiversity/data/processed/qrt_eco.parquet"
    convertParquetGeojson(file)