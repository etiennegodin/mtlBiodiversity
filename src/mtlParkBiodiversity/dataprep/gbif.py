import duckdb
from pathlib import Path 
import geopandas as gpd
import pandas as pd 
from ..core import select_file
from ..dataprep import target_crs



def ask_user_grid_file(GEOSPATIAL_PATH):

    grid_file = None
    park_file = None
    nbhood_file = None

    geospatial_files = [f for f in GEOSPATIAL_PATH.rglob("*") if f.suffix in ['.shp', '.gpkg']]

    print(f"Found {len(geospatial_files)} file in {GEOSPATIAL_PATH}")
    for file in geospatial_files:
        print(file.stem.lower())
        if "grid" in file.stem.lower():
            grid_file = file 
        elif "park" in file.stem.lower():
            park_file = file 
        elif "quartier" in file.stem.lower():
            nbhood_file = file 
    
    if grid_file is None:
        print(f'Grid file not found in {GEOSPATIAL_PATH}, please provide')
        grid_file = select_file()

    if park_file is None:
        print(f'Park file not found in {GEOSPATIAL_PATH}, please provide')
        park_file = select_file()

    if nbhood_file is None:
        print(f'Neighborhood file not found in {GEOSPATIAL_PATH}, please provide')
        nbhood_file = select_file()

    return grid_file, nbhood_file, park_file


def convert_gbif_csv(input_path, output_path, force = False, test = False, limit = None):

    if test:
        print(f'Converting to {input_path} to {output_path} file with limit set to {limit} ')
        duckdb.query(f"COPY (SELECT * FROM '{input_path}' LIMIT {limit}) TO '{output_path}' (FORMAT PARQUET)")

    else: 
        print(f'Converting to {input_path} to {output_path} file ')
        duckdb.query(f"COPY (SELECT * FROM '{input_path}') TO '{output_path}' (FORMAT PARQUET)")

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


def convert_lat_long_to_point(con, table):
    #Gbif occ to geometry 
    con.execute(f"ALTER TABLE {table} ADD COLUMN geom GEOMETRY")
    con.execute(f"UPDATE {table} SET geom = ST_Point(decimalLongitude, decimalLatitude)")

def spatial_join(gbif_occurence_db_file, park_file, output_file_path = None, test = False, limit = None, colab = False):

    # Create a connection (in-memory or persistent)
    con = duckdb.connect()  # or con = duckdb.connect("mydb.duckdb")
    if colab:
        con.execute("PRAGMA max_temp_directory_size='60GiB';")
    else:
        con.execute("PRAGMA max_temp_directory_size='25GiB';")


    # Install spatial extension 
    con.execute("INSTALL spatial;")
    con.execute("LOAD spatial;")

    #Load parks boundary file 
    con.execute(f"CREATE OR REPLACE TABLE parks AS SELECT * FROM ST_Read('{park_file}')")
    park_fields_sql = prep_park_sql(park_file)

    print('Creating gbif table...')
    #Load gbif data
    if test:
        con.execute(f"CREATE TABLE gbif AS SELECT *, ST_Point(decimalLongitude, decimalLatitude) AS geom FROM '{gbif_occurence_db_file}' LIMIT {limit}")
    else:
        if limit is not None:
            con.execute(f"CREATE TABLE gbif AS SELECT *, ST_Point(decimalLongitude, decimalLatitude) AS geom FROM '{gbif_occurence_db_file}' LIMIT {limit}")
        else:
            con.execute(f"CREATE TABLE gbif AS SELECT *, ST_Point(decimalLongitude, decimalLatitude) AS geom FROM '{gbif_occurence_db_file}'")
    
    # Add col for bbox 
    con.execute(f"ALTER TABLE gbif ADD COLUMN minx DOUBLE;")
    con.execute(f"ALTER TABLE gbif ADD COLUMN miny DOUBLE;")
    con.execute(f"ALTER TABLE gbif ADD COLUMN maxx DOUBLE;")
    con.execute(f"ALTER TABLE gbif ADD COLUMN maxy DOUBLE;")
    # Fill bbox col from geom
    con.execute(f"""UPDATE gbif 
                    SET minx = ST_XMin(geom),
                        miny = ST_YMin(geom),
                        maxx = ST_XMax(geom),
                        maxy = ST_YMax(geom);""")

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
                

                p.geom AS park_geom,
                FROM gbif g
                LEFT JOIN parks p
                    ON g.maxx >= ST_XMIN(p.geom)
                    AND g.minx <= ST_XMAX(p.geom)
                    AND g.maxy >= ST_YMIN(p.geom)
                    AND g.miny <= ST_YMAX(p.geom)
                    AND ST_Within(g.geom, p.geom) -- Spatial join predicate
                ;
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
        print(f'Successfuly saved {output_file_path}')
        con.close()
    except Exception as e:
        print(f'Failed to saved spatial join : {e}')
    return True

def prep_gbif(force = False, test = False, colab = False, limit = None):
    if colab:
        RAW_DATA_PATH = Path("/content/gdrive/MyDrive/mtlParkBiodiversity/data/raw/gbif")
        OUTPUT_PATH = Path("/content/gdrive/MyDrive/mtlParkBiodiversity/data/interim/gbif")
        GEOSPATIAL_PATH = Path("/content/mtlParkBiodiversity/data/interim/geospatial")
    else:
        RAW_DATA_PATH = Path("data/raw/gbif")
        OUTPUT_PATH = Path("data/interim/gbif")
        GEOSPATIAL_PATH = Path("data/interim/geospatial")

        Path.mkdir(OUTPUT_PATH, exist_ok= True)
        gbif_occurence_raw_file = [f for f in RAW_DATA_PATH.rglob("*.csv")][0]  # Assuming there's only one .csv file for the gbif data
    
    


    x = ask_user_grid_file(GEOSPATIAL_PATH)


    # Create output directory if not existing
    park_file =  [f for f in GEOSPATIAL_PATH.rglob("*.shp")][0]  # Assuming there's only one .shp file for the park data

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
        spatial_join(gbif_occurence_db_file,park_file, output_file_path = output_file_path, test = test, limit = limit, colab = colab)
    else:
        print("Spatial Join already done, skipping")
    
if __name__ == "__main__":

    prep_gbif(force= True)