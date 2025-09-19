import duckdb
from pathlib import Path 
import geopandas as gpd
import pandas as pd 
from ..core import select_file
from ..dataprep import target_crs

# GLOBAL VAR

# Create connection and store as class object
# Get's initialise first time it's called, otherwise object is just passed along
class DuckDBConnection:
    _instance = None

    @classmethod
    def get_connection(cls):
        if cls._instance is None:
            cls._instance = duckdb.connect("data\db\gbif.duckdb")
        return cls._instance

def check_table(table_name :str = None):

    con = DuckDBConnection.get_connection()
    df = con.execute(f"SELECT * FROM {table_name}").df()
    print(df)

def find_geospatial_fies(GEOSPATIAL_PATH: Path, debug :bool= False):

    grid_file = None
    park_file = None
    nbhood_file = None

    geospatial_files = [f for f in GEOSPATIAL_PATH.rglob("*") if f.suffix in ['.shp', '.gpkg']]

    print(f"Found {len(geospatial_files)} file in {GEOSPATIAL_PATH}")
    for file in geospatial_files:
        if debug:
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

def create_shp_file_field_sql(shp_file : Path = None , alias :str = 'p'):

    if shp_file is not None:
        sql_query = """"""

        #Get fields
        gdf = gpd.read_file(shp_file)
        columns = gdf.columns.to_list()

        #Remove geometry col by default
        if 'geometry' in columns:
            columns.remove('geometry')

        #Create string 
        for idx, col in enumerate(columns):
            x = f"""\t{alias}.{col},\n"""
            sql_query += x 
        
        return sql_query
    else:
        raise FileExistsError(f"{shp_file} is not defined")

def convert_lat_long_to_point(con, table):
    #Gbif occ to geometry 
    con.execute(f"ALTER TABLE {table} ADD COLUMN geom GEOMETRY")
    con.execute(f"UPDATE {table} SET geom = ST_Point(decimalLongitude, decimalLatitude)")
    return True

def set_geom_bbox(table_name = None):
    con = DuckDBConnection.get_connection()

    try:
        # Add col for bbox 
        con.execute(f"ALTER TABLE {table_name} ADD COLUMN minx DOUBLE;")
        con.execute(f"ALTER TABLE {table_name} ADD COLUMN miny DOUBLE;")
        con.execute(f"ALTER TABLE {table_name} ADD COLUMN maxx DOUBLE;")
        con.execute(f"ALTER TABLE {table_name} ADD COLUMN maxy DOUBLE;")
        # Fill bbox col from geom
        con.execute(f"""UPDATE {table_name} 
                            SET minx = ST_XMin(geom),
                                miny = ST_YMin(geom),
                                maxx = ST_XMax(geom),
                            maxy = ST_YMax(geom);""")
        return True
    except Exception as e:
        print(f'Could not set bbox for table {table_name}: {e}')
        return False



def create_gbif_table(gbif_occurence_db_file :Path = None, table_name = None, limit :int = None, test :bool = False):
    #Redeclare connection variable
    con = DuckDBConnection.get_connection()

    print('Creating gbif_observations table...')
    #Load gbif data

    #con.execute(f"CREATE OR REPLACE TABLE observations AS SELECT *, FROM '{gbif_occurence_db_file}'")
    if gbif_occurence_db_file is not None:
        try:
            if test and limit is not None:
                con.execute(f"CREATE OR REPLACE TABLE {table_name} AS SELECT *, ST_Point(decimalLongitude, decimalLatitude) AS geom FROM '{gbif_occurence_db_file}' LIMIT {limit}")
            else:
                if limit is not None:
                    con.execute(f"CREATE OR REPLACE TABLE {table_name} AS SELECT *, ST_Point(decimalLongitude, decimalLatitude) AS geom FROM '{gbif_occurence_db_file}' LIMIT {limit}")
                else:
                    con.execute(f"CREATE OR REPLACE TABLE {table_name} AS SELECT *, ST_Point(decimalLongitude, decimalLatitude) AS geom FROM '{gbif_occurence_db_file}'")

            set_geom_bbox(table_name= table_name)
            return True 
        
        except Exception as e:
            print('Failed to create gbif')
            print(e)
            return False
        

def create_table_from_shp(file_path : Path = None, table_name :str = None, limit : int= None, test :bool = False ):
    """
    Create duckdb tables for observations, grid, parks, nbhood
    """

    table_created = None

    #Redeclare connection variable
    con = DuckDBConnection.get_connection()
    # Install spatial extension 
    con.execute("INSTALL spatial;")
    con.execute("LOAD spatial;")

    try:
        if file_path is not None:   
            con.execute(f"CREATE OR REPLACE TABLE {table_name} AS SELECT * FROM ST_Read('{file_path}')")
            set_geom_bbox(table_name= table_name)
            table_created = True
        else:
            print('File provided to create table is None')
            table_created = False

    except Exception as e:
        print(f'Could not create table for {table_name}: {e}')
        table_created = False

    if table_created is not None:
        return table_created

def grid_spatial_join(grid_file : Path = None, output_file_path : Path = None, test : bool = False, limit :int = None, colab : bool= False):
    #Redeclare connection variable
    con = DuckDBConnection.get_connection()
    # Set final table name from expected output path name
    output_table_name = output_file_path.stem
    # Create a connection (in-memory or persistent)
    if colab:
        con.execute("PRAGMA max_temp_directory_size='60GiB';")
    else:
        con.execute("PRAGMA max_temp_directory_size='25GiB';")

    #Spatial Join 
    print("Spatial join")

    con.execute(f"""
                
                CREATE OR REPLACE TABLE {output_table_name} AS
                SELECT o.gbifID,
                o.occurrenceID,
                o.kingdom,
                o.phylum,
                o.class,
                o.order,
                o.family,
                o.genus,
                o.species,
                o.taxonRank,
                o.scientificName,
                o.eventDate,
                o.day,
                o.month,
                o.year,
                o.taxonKey,
                o.identifiedBy,
                o.basisOfRecord,
                o.license,
                o.recordedBy,
                o.geom,

                {create_shp_file_field_sql(grid_file, alias = 'g')}
                g.geom AS grid_geom,

                FROM observations o
                LEFT JOIN grid g
                    ON o.maxx >= ST_XMIN(g.geom)
                    AND o.minx <= ST_XMAX(g.geom)
                    AND o.maxy >= ST_YMIN(g.geom)
                    AND o.miny <= ST_YMAX(g.geom)
                    AND ST_Within(o.geom, g.geom) -- Spatial join predicate
                ;
                """)
    
    # Remove gbif geom point
    con.execute(f"""ALTER TABLE {output_table_name} DROP COLUMN geom;""")
    con.execute(f"""ALTER TABLE {output_table_name} RENAME COLUMN grid_geom TO geom;""")

    print('Spatial join complete, saving file...')

    con.execute('DROP TABLE IF EXISTS observations;')
    #Save 
    try:
        con.execute(f"""
                        COPY (
                            SELECT *
                            FROM {output_table_name}
                        ) TO '{output_file_path}' (FORMAT 'parquet');
                    """)
        print(f'Successfuly saved {output_file_path}')
        con.close()
        return True

    except Exception as e:
        print(f'Failed to saved spatial join : {e}')
        return False

def prep_gbif(force = False, test = False, colab = False, limit = None):

    tables_creation_dict =  {}

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
    

    grid_file, nbhood_file, park_file = find_geospatial_fies(GEOSPATIAL_PATH)
    print(grid_file, nbhood_file, park_file)
    # Create output directory if not existing
    #park_file =  [f for f in GEOSPATIAL_PATH.rglob("*.shp")][0]  # Assuming there's only one .shp file for the park data

    if test:
        print('Running gbif prep as test')
        gbif_occurence_db_file = OUTPUT_PATH / '_test_gbif_data.parquet'
        occurence_grid_file = OUTPUT_PATH / '_test_occurences_grid.parquet'

        output_file_path = OUTPUT_PATH / '_test_gbif_with_parks.parquet'
    else:
        gbif_occurence_db_file = OUTPUT_PATH / 'gbif_data.parquet'
        occurence_grid_file = OUTPUT_PATH / 'occurences_grid.parquet'

        output_file_path = OUTPUT_PATH / 'gbif_with_parks.parquet'

    # Skipping this step as parquet file is upoaded to colab directly
    if not colab:
        # Check if csv has been converted to parquet file
        if (gbif_occurence_db_file.exists() and force) or (not gbif_occurence_db_file.exists()) :
            convert_gbif_csv(gbif_occurence_raw_file, gbif_occurence_db_file, force = force, test = test, limit = limit)
        else:
            print("Convert_gbif_csv already done, skipping")

    #Iterate over files to create tables
    for file in [grid_file, nbhood_file, park_file ]:
        table_name = file.stem.split(sep= "_")[0] #Set table name as first word of file name
        created = create_table_from_shp(file_path= file, table_name=table_name , limit = limit, test = test)
        tables_creation_dict[table_name] = created

    created = create_gbif_table(gbif_occurence_db_file = gbif_occurence_db_file, table_name= 'observations', limit = limit, test = test)
    tables_creation_dict['observations'] = created

    if tables_creation_dict['grid'] and tables_creation_dict['observations']:
        if (occurence_grid_file.exists() and force) or (not occurence_grid_file.exists()):

                grid_spatial_join(grid_file = grid_file, output_file_path = occurence_grid_file,
                                test = test, limit = limit, colab = colab)
        else:
            print("Grid spatial join already done, skipping")
    
if __name__ == "__main__":

    prep_gbif(force= True)