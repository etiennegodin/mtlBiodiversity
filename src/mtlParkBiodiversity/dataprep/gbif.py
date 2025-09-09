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

output_folder = Path('data/interim/gbif/')
output_folder.mkdir(parents=True, exist_ok=True)

output_file_path = output_folder / 'gbif_with_parks.parquet'


def convert_gbif_csv(input_path, output_path, force = False, limit = None):
    if output_path.exists() and not force:
        print('File is already converted, skipping')
        return 
    
    if limit:
        print(f'Conveting to {input_path} to .parquet file with limit set to {limit} ')
        duckdb.query(f"COPY (SELECT * FROM '{input_path}' LIMIT {limit}) TO '{output_path}' (FORMAT PARQUET)")

    else: 
        print(f'Conveting to {input_path} to .parquet file ')
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

def spatial_join(gbif_occurence_db_file, park_boundaries_file, limit = 100):

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
    if limit:
        con.execute(f"CREATE TABLE gbif AS SELECT *, ST_Point(decimalLongitude, decimalLatitude) AS geom FROM '{gbif_occurence_db_file}' LIMIT {limit}")
    else:
        con.execute(f"CREATE TABLE gbif AS SELECT *, ST_Point(decimalLongitude, decimalLatitude) AS geom FROM '{gbif_occurence_db_file}'")
    
    preview_gbif_data(con, "gbif")

    #Spatial Join 
    print("Spatial join")

    con.execute(f"""
                
                CREATE OR REPLACE TABLE gbif_with_parks AS
                SELECT g.*,
                p.OBJECTID,
                p.Type,
                p.Nom,
                p.TYPO1,
                p.TYPO2
                FROM gbif g
                LEFT JOIN parks p
                ON ST_Within(g.geom, p.geom);

                """)


    #preview_gbif_data(con, "gbif_joined")
    df = con.execute("SELECT * FROM gbif_with_parks").df()
    print(df)

    print(df['OBJECTID'].unique())

    #Save 
    con.execute(f"""
    COPY gbif_with_parks TO '{output_file_path}' (FORMAT 'parquet');
    """)

    return True

def prep_gbif(force = False):

    print(f"#_{__name__}")
    # Check if csv has been converted to parquet file
    if gbif_occurence_db_file.exists and force:
        convert_gbif_csv(gbif_occurence_raw_file, gbif_occurence_db_file, force = force)
    elif not gbif_occurence_db_file.exists:
        convert_gbif_csv(gbif_occurence_raw_file, gbif_occurence_db_file, force = force)
    else:
        print("Convert_gbif_csv already done, skipping")

    # Csv to sample for tests 
    if gbif_occurence_sample_file.exists and force:
        convert_gbif_csv(gbif_occurence_raw_file, gbif_occurence_sample_file, force = force, limit = 100)
    elif not gbif_occurence_sample_file.exists:
        convert_gbif_csv(gbif_occurence_raw_file, gbif_occurence_sample_file, force = force, limit = 100)
    else:
        print("Convert_gbif_csv already done on sample, skipping")

    if output_file_path.exists and force:
        spatial_join(gbif_occurence_db_file,park_boundaries_file,limit = None)
    elif not output_file_path.exists:
        spatial_join(gbif_occurence_db_file,park_boundaries_file,limit = None)
    else:
        print("Spatial Join already done, skipping")
    
if __name__ == "__main__":
    prep_gbif(force= True)