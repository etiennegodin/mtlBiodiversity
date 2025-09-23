from pathlib import Path 
import geopandas as gpd
from mtlBio.core import read_sql_template, DuckDBConnection

# GLOBAL VAR
tables_creation_dict = {}

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
    



def grid_spatial_join2(left_table_name : str = None, right_table_name: str = None, output_file_path : Path = None, test : bool = False, limit :int = None):
    #Redeclare connection variable
    con = DuckDBConnection.get_connection()
    # Set final table name from expected output path name
    output_table_name = output_file_path.stem
    
    if check_table_exists(left_table_name) and check_table_exists(right_table_name):
        
        # Set limit for tmp files on disk
        con.execute("PRAGMA max_temp_directory_size='25GiB';")
        
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
        
        print(f"Spatial join of {left_table_name} & {right_table_name}...")

        # Main spatial join query 
        con.execute(sjoin_sql_query)
        
        clean_joined_table_template = read_sql_template('clean_gbif_sjoin')
        clean_joined_table_query = clean_joined_table_template.render(output_table_name = output_table_name)
        print('Cleaning joined table')
        try:
            con.execute(clean_joined_table_query)
        except Exception as e:
            print(e)
            
        print('Spatial join complete, saving file...')

        #Save 
        try:
            con.execute(f"""
                            COPY (
                                SELECT *
                                FROM {output_table_name}
                            ) TO '{output_file_path}' (FORMAT 'parquet');
                        """)
            print(f'Successfuly saved {output_file_path}')
            
            #Apend to dict of created tables to keep track 
            global tables_creation_dict
            tables_creation_dict[output_table_name] = True
            return output_table_name

        except Exception as e:
            print(f'Failed to saved spatial join : {e}')
            return None
        
        
    elif not check_table_exists(left_table_name) and check_table_exists(right_table_name):
        print(f"Left table {left_table_name} doesn't exists ")
        return None
    elif check_table_exists(left_table_name) and not check_table_exists(right_table_name):
        print(f"Right table {right_table_name} doesn't exists ")
        return None
    else:
        print(f"Tables {left_table_name}, {right_table_name} do not exist")
        return None



def gbif_spatial_joins(gbif_occurence_db_file :Path = None, force = False, test = False, limit = None):
    """Performs main geospatial joins of gbif observation to set of shp files

    Args:
        gbif_occurence_db_file (Path, optional): _description_. Defaults to None.
        force (bool, optional): _description_. Defaults to False.
        test (bool, optional): _description_. Defaults to False.
        limit (_type_, optional): _description_. Defaults to None.
    """
    global tables_creation_dict 

    OUTPUT_PATH = Path("data/interim/joined")
    GEOSPATIAL_PATH = Path("data/interim/geospatial")

    # Find all expected geospatial files
    #grid_file, nbhood_file, park_file = find_geospatial_fies(GEOSPATIAL_PATH)
    files = find_files(folder_path= GEOSPATIAL_PATH, expected_files= ['grid', 'quartiers', 'park'], suffix= '.shp')
    
    print(files)

    if test:
        print('Running gbif prep as test')
        occurence_grid_file = OUTPUT_PATH / '_test_occurences_grid.parquet'
        quartier_joined_file = OUTPUT_PATH / '_test_occurences_quartiers.parquet'
    else:
        occurence_grid_file = OUTPUT_PATH / 'occurences_grid.parquet'
        quartier_joined_file = OUTPUT_PATH / 'occurences_quartiers.parquet'

    #Iterate over files to create tables
    for file in files:
        table_name = file.stem.split(sep= "_")[0] #Set table name as first word of file name
        created = create_table_from_shp(file_path= file, table_name=table_name , limit = limit, test = test)
        tables_creation_dict[table_name] = created
    print(gbif_occurence_db_file)
    

    if gbif_occurence_db_file is not None:
        created = create_gbif_table(gbif_occurence_db_file = gbif_occurence_db_file, table_name= 'observations', limit = limit, test = test)
        tables_creation_dict['observations'] = created

    if tables_creation_dict['grid'] and tables_creation_dict['observations']:
        if (occurence_grid_file.exists() and force) or (not occurence_grid_file.exists()):

            grid_joined_table = grid_spatial_join2(left_table_name = 'observations', right_table_name= 'grid', output_file_path = occurence_grid_file, test  = False, limit = None)
            if grid_joined_table is not None:
                grid_spatial_join2(left_table_name = grid_joined_table, right_table_name= 'quartiers', output_file_path = quartier_joined_file, test  = False, limit = None)
            
        else:
            print("Grid spatial join already done, skipping")
    
if __name__ == "__main__":

    pass