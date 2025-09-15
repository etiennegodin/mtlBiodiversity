import duckdb
from pathlib import Path 
import geopandas as gpd
import pandas as pd 
from ..core import clip_to_region, convert_crs
from ..dataprep import target_crs





def convert_gbif_csv(input_path, output_path, force = False, test = False, limit = None):

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

def prep_park_sql(park_file, alias = 'p'):

    sql_query = """"""

    #Get fields
    gdf = gpd.read_file(park_file)
    columns = gdf.columns.to_list()

    #Remove geometry col by default
    if 'geometry' in columns:
        columns.remove('geometry')

    #Create string 
    for idx, col in enumerate(columns):
        x = f"""\t{alias}.{col},\n"""
        sql_query += x 

    return sql_query

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

def spatial_join(gbif_occurence_db_file, park_file, output_file_path = None, test = False, limit = None):

    # Create a connection (in-memory or persistent)
    con = duckdb.connect()  # or con = duckdb.connect("mydb.duckdb")

    con.execute("PRAGMA max_temp_directory_size='25GiB';")

    # Install spatial extension 
    con.execute("INSTALL spatial;")
    con.execute("LOAD spatial;")

    #Load parks boundary file 
    con.execute(f"CREATE OR REPLACE TABLE parks AS SELECT * FROM ST_Read('{park_file}')")
    park_fields_sql = prep_park_sql(park_file)

    #preview_gbif_data(con, "parks")
    print('Creating gbif table...')
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

                {park_fields_sql}

                p.geom AS park_geom
                FROM gbif g
                LEFT JOIN parks p
                ON ST_Within(g.geom, p.geom);

                """)

    print('Spatial join complete, saving file...')

    con.execute('DROP TABLE IF EXISTS gbif;')
    #Save 
    try:
        con.execute(f"""
                        COPY (
                            SELECT *
                            FROM gbif_with_parks
                        ) TO '{output_file_path}' (FORMAT 'parquet');
                    """)
        con.close()
    except Exception as e:
        print(f'Failed to saved spatial join : {e}')
    return True

def prep_gbif(force = False, test = False, colab = False, limit = None):
    if colab:
        RAW_DATA_PATH = Path("/content/gdrive/MyDrive/mtlParkBiodiversity/data/raw/gbif")
        OUTPUT_PATH = Path("/content/gdrive/MyDrive/mtlParkBiodiversity/data/interim/gbif")
        PARK_PATH = Path("/content/gdrive/MyDrive/mtlParkBiodiversity/data/interim/parks")
    else:
        RAW_DATA_PATH = Path("data/raw/gbif")
        OUTPUT_PATH = Path("data/interim/gbif")
        PARK_PATH = Path("data/interim/parks")

        Path.mkdir(OUTPUT_PATH, exist_ok= True)
        gbif_occurence_raw_file = [f for f in RAW_DATA_PATH.rglob("*.csv")][0]  # Assuming there's only one .csv file for the gbif data


    # Create output directory if not existing
    park_file =  [f for f in PARK_PATH.rglob("*.shp")][0]  # Assuming there's only one .shp file for the park data

    if test:
        print('Running gbif prep as test')
        gbif_occurence_db_file = OUTPUT_PATH / '_test_gbif_data.parquet'
        output_file_path = OUTPUT_PATH / '_test_gbif_with_parks.parquet'
    else:
        gbif_occurence_db_file = OUTPUT_PATH / 'gbif_data.parquet'
        output_file_path = OUTPUT_PATH / 'gbif_with_parks.parquet'

    # Skipping this step as parquet file is upoaded to colab directly
    if not colab:
        # Check if csv has been converted to parquet file
        if (gbif_occurence_db_file.exists() and force) or (not gbif_occurence_db_file.exists()) :
            convert_gbif_csv(gbif_occurence_raw_file, gbif_occurence_db_file, force = force, test = test, limit = limit)
        else:
            print("Convert_gbif_csv already done, skipping")

    if (output_file_path.exists() and force) or (not output_file_path.exists()):
        pass
        spatial_join(gbif_occurence_db_file,park_file, output_file_path = output_file_path, test = test, limit = limit)
    else:
        print("Spatial Join already done, skipping")
    
if __name__ == "__main__":
    prep_gbif(force= True)