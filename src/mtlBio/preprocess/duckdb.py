from pathlib import Path 
from mtlBio.core import DuckDBConnection, convertToPath, assign_table_alias, load_spatial_extension, set_geom_bbox



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

    except Exception as e:
        print(f'Could not create table for {table_name}: {e}')
    
def create_gbif_table(gbif_data_path :Path = None,  marker_file:str = None):

    db = DuckDBConnection()
    con = db.conn
    
    # Install spatial extension 
    load_spatial_extension(con)
    
    gbif_data_path = convertToPath(gbif_data_path)
    table_name = 'gbif_raw'
        
    print('Creating gbif_observations table...')

    if gbif_data_path is not None:
        
        query = f"""CREATE OR REPLACE TABLE {table_name} AS
                SELECT *,
                ST_Point(decimalLongitude, decimalLatitude) AS geom,
                FROM '{gbif_data_path}'
                """
        
        #print(query)
        try:
            con.execute(query)
            set_geom_bbox(table_name= table_name)
            Path(marker_file).touch()
            
        except Exception as e:
            print('Failed to create gbif')
            print(e)
    

def filter_gbif_data(gbif_cols:str = None, limit :int = None,  marker_file:str = None,  coordUncerFilter:float = None, wkt_filters: str = None ):
    
    table_name = 'filtered_gbif'
    
    if limit == 'None':
        limit = None
    
    db = DuckDBConnection()
    con = db.conn
    
    # Install spatial extension 
    load_spatial_extension(con)
    
    #COnvert cols if not list 
    if not isinstance(gbif_cols, list):
        gbif_cols_temp = gbif_cols.split(sep= ",")
        gbif_cols = []
        for col in gbif_cols_temp:
            col = col.strip()
            gbif_cols.append(col)
     
    gbif_cols = assign_table_alias(gbif_cols, alias= 'g')
             
    query = f"""CREATE OR REPLACE TABLE {table_name} AS
                SELECT {gbif_cols}
                FROM gbif_raw AS g,
                WHERE coordinateUncertaintyInMeters <= {coordUncerFilter}
                {'LIMIT ' + str(limit) if (limit is not None) else ''}
                """ 
                
    #print(query)
    try:
        con.execute(query)
        set_geom_bbox(table_name= table_name)
        Path(marker_file).touch()
            
    except Exception as e:
        print('Failed to create gbif')
        print(e)
    
