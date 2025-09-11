import duckdb
from pathlib import Path 
import geopandas as gpd
import pandas as pd 
from ..core import clip_to_region, convert_crs
from ..dataprep import target_crs

RAW_DATA_PATH = Path("data/raw/gbif")
OUTPUT_PATH = Path("data/interim/gbif")
PARK_FILE_PATH = Path("data/interim/Espace_Vert_clipped.shp")

def convert_gbif_csv(input_path, output_path, force = False, test = False, limit = None):
    if output_path.exists() and not force:
        print('File is already converted, skipping')
        return 
    
    if test:
        print(f'Converting to {input_path} to .parquet file with limit set to {limit} ')
        duckdb.query(f"COPY (SELECT * FROM '{input_path}' LIMIT {limit}) TO '{output_path}' (FORMAT PARQUET)")

    else: 
        print(f'Converting to {input_path} to .parquet file ')
        duckdb.query(f"COPY (SELECT * FROM '{input_path}') TO '{output_path}' (FORMAT PARQUET)")
  

def preview_gbif_data(con, table, limit = 100):
    print('-'*25, 'Preview', '-'*25)
    df = con.execute(f"SELECT * FROM {table} LIMIT {limit}").df()
    print(df.head(limit))
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

def spatial_join(gbif_occurence_db_file, PARK_FILE_PATH, test = False, limit = None):

    # Create a connection (in-memory or persistent)
    con = duckdb.connect()  # or con = duckdb.connect("mydb.duckdb")

    con.execute("PRAGMA max_temp_directory_size='25GiB';")



    # Install spatial extension 
    con.execute("INSTALL spatial;")
    con.execute("LOAD spatial;")

    #Load parks boundary file 
    con.execute(f"CREATE OR REPLACE TABLE parks AS SELECT * FROM ST_Read('{PARK_FILE_PATH}')")

    #preview_gbif_data(con, "parks")

    #Load gbif data
    if test:
        con.execute(f"CREATE TABLE gbif AS SELECT *, ST_Point(decimalLongitude, decimalLatitude) AS geom FROM '{gbif_occurence_db_file}' LIMIT {limit}")
    else:
        con.execute(f"CREATE TABLE gbif AS SELECT *, ST_Point(decimalLongitude, decimalLatitude) AS geom FROM '{gbif_occurence_db_file}'")
    
    #preview_gbif_data(con, "gbif")

    #Spatial Join 
    print("Spatial join")

    con.execute(f"""
                
                CREATE OR REPLACE TABLE gbif_with_parks AS
                SELECT g.gbifID,
                g.occurrenceID,
                g.kingdom,
                g.phylum,
                g.class,
                g.order,
                g.family,
                g.genus,
                g.species,
                g.taxonRank,
                g.scientificName,
                g.eventDate,
                g.day,
                g.month,
                g.year,
                g.taxonKey,
                g.identifiedBy,
                g.basisOfRecord,
                g.license,
                g.recordedBy,
                g.geom,

                # create dynamic query 
                (SELECT * EXCEPT xxx) is not possible 
                # read park.shp and get columns, remove geom, (p.park_id, p.park_name, etc) - have it dynamic without p.geom to avoid duplicates 
                # add , at end

                # handle geom 
                p.geom AS park_geom
                FROM gbif g
                LEFT JOIN parks p
                ON ST_Within(g.geom, p.geom);

                """)

    print('Spatial join complete, saving file...')

    con.execute('DROP TABLE IF EXISTS gbif;')
    con.execute('PRAGMA optimize;')
    print('Dropped gbif table to save memory')
    #Save 
    #con.execute(f"""COPY gbif_with_parks TO '{output_file_path}' (FORMAT 'parquet');""")
    con.execute(f"""
                    COPY (
                        SELECT *
                        FROM gbif_with_parks
                    ) TO '{output_file_path}' (FORMAT 'parquet');
                """)
    con.close()
    return True

def prep_gbif(force = False, test = False, limit = None):

    # Create output directory if not existing
    Path.mkdir(OUTPUT_PATH, exist_ok= True)

    gbif_occurence_raw_file = [f for f in RAW_DATA_PATH.rglob("*.csv")][0]  # Assuming there's only one .csv file for the gbif data

    if test:
        gbif_occurence_db_file = OUTPUT_PATH / '_gbif_data_test.parquet'
        output_file_path = OUTPUT_PATH / '_gbif_with_parks_test.parquet'
    else:
        gbif_occurence_db_file = OUTPUT_PATH / 'gbif_data.parquet'
        output_file_path = OUTPUT_PATH / 'gbif_with_parks.parquet'


    print(f"#_{__name__}")
    # Check if csv has been converted to parquet file
    if (gbif_occurence_db_file.exists and force) or (not gbif_occurence_db_file.exists) :
        convert_gbif_csv(gbif_occurence_raw_file, gbif_occurence_db_file, force = force, test = test, limit = limit)
    else:
        print("Convert_gbif_csv already done, skipping")

    if (output_file_path.exists and force) or (not output_file_path.exists):
        spatial_join(gbif_occurence_db_file,PARK_FILE_PATH, test = test, limit = limit)
    else:
        print("Spatial Join already done, skipping")
    
if __name__ == "__main__":
    prep_gbif(force= True)