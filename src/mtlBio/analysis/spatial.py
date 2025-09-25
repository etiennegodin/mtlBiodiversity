from pathlib import Path 
from mtlBio.core import read_sql_template, DuckDBConnection, assign_table_alias, load_spatial_extension


def get_table_columns(table_name = None, con = None):
    columns = None
    if table_name is not None:
        if con is not None:
            columns = con.execute(f"PRAGMA table_info('{table_name}')").df()
            columns = columns['name'].tolist()
        
    if isinstance(columns, list):
        return columns


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

    
