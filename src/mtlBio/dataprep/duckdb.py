from pathlib import Path 
import geopandas as gpd
from mtlBio.core import read_sql_template, DuckDBConnection, convertToPath

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




def create_table_from_shp(file_path : Path = None, limit : int= None):
    """
    Create duckdb tables for observations, grid, parks, nbhood
    """
    #Redeclare connection variable
    con = DuckDBConnection.get_connection()
    # Install spatial extension 
    con.execute("INSTALL spatial;")
    con.execute("LOAD spatial;")
    
    file_path = convertToPath(file_path)
    
    #Define table name from file
    table_name = file_path.stem.split(sep= '_')[0]

    try:
        if file_path is not None:   
            con.execute(f"CREATE OR REPLACE TABLE {table_name} AS SELECT * FROM ST_Read('{file_path}')")
            set_geom_bbox(table_name= table_name)
        else:
            print('File provided to create table is None')

    except Exception as e:
        print(f'Could not create table for {table_name}: {e}')


def create_gbif_table(gbif_data_path :Path = None, limit :int = None):
    
    gbif_data_path = convertToPath(gbif_data_path)
    table_name = gbif_data_path.stem
    
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
                """
    
    #Redeclare connection variable
    con = DuckDBConnection.get_connection()
    
    print('Creating gbif_observations table...')
    #Load gbif data

    if gbif_data_path is not None:
        
        query = f"""CREATE OR REPLACE TABLE {table_name} AS
                SELECT {gbif_raw_col}
                ST_Point(decimalLongitude, decimalLatitude) AS geom,
                FROM '{gbif_data_path}' AS f
                {'LIMIT ' + str(limit) if (limit is not None) else ''} """
        
        try:
            con.execute(query)
            set_geom_bbox(table_name= table_name)
            return True 
        
        except Exception as e:
            print('Failed to create gbif')
            print(e)
            return False
        