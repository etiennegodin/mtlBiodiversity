import duckdb
from pathlib import Path 
import geopandas as gpd
import pandas as pd 
from ..core import clip_to_region, convert_crs
from ..dataprep import target_crs

gbif_folder = Path('C:/Users/manat/Documents/Projects/mtlParkBiodiversity/data/raw/gbif/')

gbif_occurence_raw_file = gbif_folder / 'gbif_occurences.csv'
gbif_occurence_db_file = gbif_folder / 'gbif_occurences.parquet'
gbif_occurence_sample_file = gbif_folder / 'gbif_occurences_sample.parquet'

park_boundaries_file = Path("data/interim/Espace_Vert_clipped.shp")


def convert_gbif_csv(input_path, output_path, force = False, limit = None):
    if output_path.exists() and not force:
        print('File is already converted, skipping')
        return 
    
    if limit:
        duckdb.query(f"COPY (SELECT * FROM '{input_path}' LIMIT {limit}) TO '{output_path}' (FORMAT PARQUET)")
    else: 
        duckdb.query(f"COPY (SELECT * FROM '{input_path}') TO '{output_path}' (FORMAT PARQUET)")
  

def preview_gbif_data(con, table, limit = 100):
    print('-'*25, 'Preview', '-'*25)
    df = con.execute(f"SELECT * FROM {table} LIMIT {limit}").df()
    print(df.head())
    print(df.columns)
    print('/'*50)


def check_crs(file, debug = False):
    if debug:
        print(file.suffix)
    if file.suffix == ".shp":
        try:
            x = gpd.read_file(file)
            crs = x.crs
            print(crs)
            return crs

        except Exception as e:
            print(f'Error reading {file} : {e}')
            return False
    elif file.suffix == '.parquet':
        try:
            crs = x.crs
            print(crs)
            return crs

        except Exception as e:
            print(f'Error reading {file} : {e}')
            return False
    else:
        print(f'Error reading {file}')
        return False
    
    
def get_gbif_data_crs(gbif_occurence_raw_file, lat, long, limit = 10, ):
    from shapely.geometry import Point

    df = pd.read_csv(gbif_occurence_raw_file, nrows = limit, delimiter = '\t')
    gdf = gpd.GeoDataFrame(df, geometry=[Point(xy) for xy in zip(df["decimalLongitude"], df["decimalLatitude"])] , crs = 4326 )
    print(gdf.head())
    print(gdf.crs)
    pass

def convert_lat_long_to_point(con, table):
    #Gbif occ to geometry 
    con.execute(f"ALTER TABLE {table} ADD COLUMN geom GEOMETRY")
    con.execute(f"UPDATE {table} SET geom = ST_Point(decimalLongitude, decimalLatitude)")

def spatial_join(limit = 100):

    # Create a connection (in-memory or persistent)
    con = duckdb.connect()  # or con = duckdb.connect("mydb.duckdb")

    # Install spatial extension 
    con.execute("INSTALL spatial;")
    con.execute("LOAD spatial;")

    #Load parks boundary file 
    con.execute(f"CREATE OR REPLACE TABLE parks AS SELECT * FROM ST_Read('{park_boundaries_file}')")
    print(con.execute("DESCRIBE parks").fetchall())

    preview_gbif_data(con, "parks")

    #Load gbif data
    con.execute(f"CREATE TABLE gbif AS SELECT *, ST_Point(decimalLongitude, decimalLatitude) AS geom FROM '{gbif_occurence_db_file}' LIMIT {limit}")

    preview_gbif_data(con, "gbif")

    print("Spatial join")
    # Create new col for park id 
    #Spatial Join 

    con.execute(f"""CREATE OR REPLACE TABLE gbif_with_parks AS
                SELECT g.*, p.OBJECTID
                FROM gbif g
                LEFT JOIN parks p
                ON ST_Within(g.geom, p.geom);
                """)


    #preview_gbif_data(con, "gbif_joined")
    df = con.execute("SELECT * FROM gbif_with_parks").df()
    print(df)

    print(df['OBJECTID'].unique())
    

    print('d')

    

    return

    con.execute(f"""
                    
                    
                    ALTER TABLE gbif ADD COLUMN Nom VARCHAR;
                    UPDATE gbif
                    SET Nom = (
                    SELECT p.Nom
                    FROM parks p
                    WHERE ST_Within(p.geom,gbif.geom )
                    );
                    
                    
                    """)

def prep_gbif(force = False):
    print(f"#_{__name__}")
    convert_gbif_csv(gbif_occurence_raw_file, gbif_occurence_db_file, force = force)
    convert_gbif_csv(gbif_occurence_raw_file, gbif_occurence_sample_file, force = force, limit = 100)

    #get_gbif_data_crs(gbif_occurence_raw_file, 'decimalLatitude', 'decimalLongitude')
    print('Check crs')
    check_crs(park_boundaries_file, debug = False)

    spatial_join(limit = 10)

if __name__ == "__main__":
    prep_gbif(force= True)