from pathlib import Path 
from mtlBio.core import read_sql_template, DuckDBConnection, convertToPath

def load_spatial_extension(con = None):
    if con is None:
        db = DuckDBConnection()
        con = db.conn
        
    con.execute("INSTALL spatial;")
    con.execute("LOAD spatial;")

def assign_table_alias(columns: list = None, alias :str = None):
    query = """"""
    for col in columns:
        col_new = f'{alias}.{col}'
        query += f"{col_new},\n\t\t\t\t"
    return query

def get_table_columns(table_name = None, con = None):
    columns = None
    if table_name is not None:
        if con is not None:
            columns = con.execute(f"PRAGMA table_info('{table_name}')").df()
            columns = columns['name'].tolist()
        
    if isinstance(columns, list):
        return columns

def set_geom_bbox(table_name = None):
    db = DuckDBConnection()
    con = db.conn
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



def create_table_from_shp(file_path : Path = None, marker_file:str = None):
    """
    Create duckdb tables for observations, grid, parks, nbhood
    """
    #Redeclare connection variable
    db = DuckDBConnection()
    con = db.conn
    
    # Install spatial extension 
    load_spatial_extension(con)
    
    file_path = convertToPath(file_path)
    
    #Define table name from file
    table_name = file_path.stem.split(sep= '_')[0]

    try:
        if file_path is not None:   
            con.execute(f"CREATE OR REPLACE TABLE {table_name} AS SELECT * FROM ST_Read('{file_path}')")
            set_geom_bbox(table_name= table_name)
            Path(marker_file).touch()

            print('File provided to create table is None')

    except Exception as e:
        print(f'Could not create table for {table_name}: {e}')
    
def create_gbif_table(gbif_data_path :Path = None, limit :int = None, marker_file:str = None, coordUncerFilter:float = None ):
    
    if limit == 'None':
        limit = None
    
    db = DuckDBConnection()
    con = db.conn
    
    # Install spatial extension 
    load_spatial_extension(con)

    gbif_data_path = convertToPath(gbif_data_path)
    table_name = gbif_data_path.stem.split(sep='_')[0]
    
    gbif_raw_col = """f.gbifID,
                f.occurrenceID,
                f.kingdom,
                f.phylum,
                f.class,
                f.order,
                f.family,
                f.genus,
                f.species,
                f.taxonRank,
                f.scientificName,
                f.eventDate,
                f.day,
                f.month,
                f.year,
                f.taxonKey,
                f.basisOfRecord,
                f.license,
                f.recordedBy,
                f.issue,
                f.publishingOrgKey,
                f.coordinateUncertaintyInMeters,
                """
    
    #Redeclare connection variable
    
    print('Creating gbif_observations table...')
    #Load gbif data

    if gbif_data_path is not None:
        
        query = f"""CREATE OR REPLACE TABLE {table_name} AS
                SELECT {gbif_raw_col}
                ST_Point(decimalLongitude, decimalLatitude) AS geom,
                FROM '{gbif_data_path}' AS f
                WHERE coordinateUncertaintyInMeters <= {coordUncerFilter}
                {'LIMIT ' + str(limit) if (limit is not None) else ''} """
        
        try:
            con.execute(query)
            set_geom_bbox(table_name= table_name)
            Path(marker_file).touch()
            
        
        except Exception as e:
            print('Failed to create gbif')
            print(e)
    

    

def grid_spatial_join(left_table_name : str = None, right_table_name: str = None, marker_file:str= None):
    #Redeclare connection variable
    db = DuckDBConnection()
    con = db.conn

    load_spatial_extension(con = con )
    # Set final table name from expected output path name
    output_table_name = f"{right_table_name}_sjoin"
    
        
    # Set limit for tmp files on disk
    con.execute("PRAGMA max_temp_directory_size='25GiB';")
    
    left_table_cols = get_table_columns(table_name= left_table_name, con = con)
    right_table_cols= get_table_columns(table_name= right_table_name, con = con)
    
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
    
    print(f"Spatial join of {left_table_name} & {right_table_name}...")

    # Main spatial join query 
    try:
        con.execute(sjoin_sql_query)
        
        clean_joined_table_template = read_sql_template('clean_gbif_sjoin')
        clean_joined_table_query = clean_joined_table_template.render(output_table_name = output_table_name)
    
        print('Cleaning joined table')
        try:
            con.execute(clean_joined_table_query)
            print('Spatial join complete')
            Path(marker_file).touch()

        except Exception as e:
            print(e)
    except Exception as e:
        print(e)    

    
