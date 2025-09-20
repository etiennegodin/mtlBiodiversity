import duckdb
from pathlib import Path 
import geopandas as gpd
import pandas as pd 
from mtlBio.core import select_file, read_sql_template
from mtlBio.dataprep import target_crs
from jinja2 import Template

# GLOBAL VAR

# Create connection and store as class object
# Get's initialise first time it's called, otherwise object is just passed along
class DuckDBConnection:
    _instance = None

    @classmethod
    def get_connection(cls):
        if cls._instance is None:
            cls._instance = duckdb.connect("data\db\gbif_spatial_join.duckdb")
        return cls._instance


def check_table_exists(table_name = None):
    con = DuckDBConnection.get_connection()
  
    if table_name is not None:
        result = con.execute(f"""
            SELECT COUNT(*) 
            FROM information_schema.tables 
            WHERE table_name = '{table_name.lower()}'
        """).fetchone()
    if result[0] > 0:
        return True
    else:
        return False
    
def inspect_table(table_name :str = None):

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

    query = None
    
    if test:
        print(f'Converting to {input_path} to {output_path} file with limit set to {limit} ')
        query = f"COPY (SELECT * FROM '{input_path}' LIMIT {limit}) TO '{output_path}' (FORMAT PARQUET)"
    else: 
        print(f'Converting to {input_path} to {output_path} file ')
        query = f"COPY (SELECT * FROM '{input_path}') TO '{output_path}' (FORMAT PARQUET)"
        
    if query is not None:
        try:
            duckdb.query(query)
            return True
        except Exception as e:
            print(f"Failed to convert csv to parquet: {e}")
            return False
        
def assign_table_alias(columns: list = None, alias :str = None):
    query = """"""
    for col in columns:
        col_new = f'{alias}.{col}'
        query += f"{col_new},\n\t\t\t\t"
    return query
        

def get_table_columns(table_name = None):
    columns = None
    if table_name is not None:
        con = DuckDBConnection.get_connection()
        columns = con.execute(f"PRAGMA table_info('{table_name}')").df()
        columns = columns['name'].tolist()
        
    if isinstance(columns, list):
        return columns
    
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
    
    
def grid_spatial_join2(left_table_name : str = None, right_table_name: str = None, output_file_path : Path = None, test : bool = False, limit :int = None):
    #Redeclare connection variable
    con = DuckDBConnection.get_connection()
    # Set final table name from expected output path name
    output_table_name = output_file_path.stem
    
    if check_table_exists(left_table_name) and check_table_exists(right_table_name):
        
        # SEt limit for tmp files on disk
        con.execute("PRAGMA max_temp_directory_size='25GiB';")

        #Spatial Join 
        print("Grid spatial join")
        
        left_table_cols = get_table_columns(table_name= left_table_name)
        right_table_cols= get_table_columns(table_name= right_table_name)
        
        # Remove right col as it's handled by the sjoin sql query ( renamed to right geom temporrily to avoid confusion)
        if 'geom' in right_table_cols:
            right_table_cols.remove('geom')
        
        left_table_cols_alias = assign_table_alias(left_table_cols, alias = 'l')
        right_table_cols_alias = assign_table_alias(right_table_cols, alias = 'r')
       
       
        sjoin_sql_template= read_sql_template('gbif_spatial_join')
        sjoin_sql_query = sjoin_sql_template.render(output_table_name =output_table_name,
                               left_fields = left_table_cols_alias,
                               right_fields = right_table_cols_alias,
                               left_table_name = left_table_name,
                               right_table_name = right_table_name

                               )
       
        with open('mysql.sql', 'w') as f:
           f.write(sjoin_sql_query)
       
       
        con.execute(sjoin_sql_query)
        
        # Remove gbif geom point
        con.execute(f"""ALTER TABLE {output_table_name} DROP COLUMN geom;""")
        con.execute(f"""ALTER TABLE {output_table_name} RENAME COLUMN right_geom TO geom;""")

        print('Spatial join complete, saving file...')

        #con.execute('DROP TABLE IF EXISTS observations;')
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
        
        
    elif not check_table_exists(left_table_name) and check_table_exists(right_table_name):
        print(f"Left table {left_table_name} doesn't exists ")
        return False
    elif check_table_exists(left_table_name) and not check_table_exists(right_table_name):
        print(f"Right table {right_table_name} doesn't exists ")
        return False
    else:
        print(f"Tables {left_table_name}, {right_table_name} do not exist")
        return False




def grid_spatial_join(grid_file : Path = None, output_file_path : Path = None, test : bool = False, limit :int = None):
    #Redeclare connection variable
    con = DuckDBConnection.get_connection()
    # Set final table name from expected output path name
    output_table_name = output_file_path.stem

    # SEt limit for tmp files on disk
    con.execute("PRAGMA max_temp_directory_size='25GiB';")

    #Spatial Join 
    print("Grid spatial join")

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


def prep_gbif_data(force = False, test = False, limit = None):
    """Create .parquet file from gbif obsrvations csv file

    Args:
        force (bool, optional): _description_. Defaults to False.
        test (bool, optional): _description_. Defaults to False.
        limit (_type_, optional): _description_. Defaults to None.
    """
    RAW_DATA_PATH = Path("data/raw/gbif")
    OUTPUT_PATH = Path("data/interim/gbif")
    gbif_occurence_db_file = None

    if test:
        print('Running gbif prep as test')
        gbif_occurence_db_file = OUTPUT_PATH / '_test_gbif_data.parquet'
    else:
        gbif_occurence_db_file = OUTPUT_PATH / 'gbif_data.parquet'

    gbif_occurence_raw_file = [f for f in RAW_DATA_PATH.rglob("*.csv")][0]  # Assuming there's only one .csv file for the gbif data

    # Check if csv has been converted to parquet file
    if (gbif_occurence_db_file.exists() and force) or (not gbif_occurence_db_file.exists()) :
        converted = convert_gbif_csv(gbif_occurence_raw_file, gbif_occurence_db_file, force = force, test = test, limit = limit)
        if converted:
            return gbif_occurence_db_file
        else:
            return None
    else:
        print("Convert_gbif_csv already done, skipping")
        return gbif_occurence_db_file

def gbif_spatial_joins(gbif_occurence_db_file :Path = None, force = False, test = False, limit = None):
    """Performs main geospatial joins of gbif observation to set of shp files

    Args:
        gbif_occurence_db_file (Path, optional): _description_. Defaults to None.
        force (bool, optional): _description_. Defaults to False.
        test (bool, optional): _description_. Defaults to False.
        limit (_type_, optional): _description_. Defaults to None.
    """
    tables_creation_dict =  {}

    OUTPUT_PATH = Path("data/interim/joined")
    GEOSPATIAL_PATH = Path("data/interim/geospatial")

    # Find all expected geospatial files
    grid_file, nbhood_file, park_file = find_geospatial_fies(GEOSPATIAL_PATH)
    print(grid_file, nbhood_file, park_file)

    if test:
        print('Running gbif prep as test')
        occurence_grid_file = OUTPUT_PATH / '_test_occurences_grid.parquet'
        quartier_joined_file = OUTPUT_PATH / '_test_gbif_with_parks.parquet'
    else:
        occurence_grid_file = OUTPUT_PATH / 'occurences_grid.parquet'
        quartier_joined_file = OUTPUT_PATH / 'gbif_with_parks.parquet'

    #Iterate over files to create tables
    for file in [grid_file, nbhood_file, park_file ]:
        table_name = file.stem.split(sep= "_")[0] #Set table name as first word of file name
        created = create_table_from_shp(file_path= file, table_name=table_name , limit = limit, test = test)
        tables_creation_dict[table_name] = created
    print(gbif_occurence_db_file)
    

    if gbif_occurence_db_file is not None:
        created = create_gbif_table(gbif_occurence_db_file = gbif_occurence_db_file, table_name= 'observations', limit = limit, test = test)
        tables_creation_dict['observations'] = created

    if tables_creation_dict['grid'] and tables_creation_dict['observations']:
        if (occurence_grid_file.exists() and force) or (not occurence_grid_file.exists()):

                grid_spatial_join2(left_table_name = 'observations', right_table_name= 'grid', output_file_path = occurence_grid_file, test  = False, limit = None)
        else:
            print("Grid spatial join already done, skipping")
    
if __name__ == "__main__":

    pass