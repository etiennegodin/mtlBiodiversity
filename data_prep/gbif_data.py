import duckdb
from pathlib import Path 
gbif_occurence_raw_file = Path('C:/Users/manat/Documents/Projects/mtlParkBiodiversity/data/raw/gbif/gbif_occurences.csv')
gbif_occurence_db_file = Path('C:/Users/manat/Documents/Projects/mtlParkBiodiversity/data/raw/gbif/gbif_occurences.parquet')

park_boundaries_file = "data/raw/espace_verts/Espace_Vert.shp"

def convert_gbif_csv(input_path, output_path, overwrite = False):
    if output_path.exists() and not overwrite:
        print('File is already converted, skipping')
        return 
    
    duckdb.query(f"COPY (SELECT * FROM '{input_path}') TO '{output_path}' (FORMAT PARQUET)")

def preview_gbif_data(con, table, limit = 100):
    print('-'*25, 'Preview', '-'*25)
    df = con.execute(f"SELECT * FROM {table} LIMIT {limit}").df()
    print(df.head())
    print(df.columns)
    print()




def spatial_join():

    # Create a connection (in-memory or persistent)
    con = duckdb.connect()  # or con = duckdb.connect("mydb.duckdb")

    # Install spatial extension 
    con.execute("INSTALL spatial; LOAD spatial;")

    #Load parks boundary file 
    con.execute(f"CREATE TABLE parks AS SELECT * FROM ST_Read('{park_boundaries_file}')")

    #Load gbif data
    con.execute(f"CREATE TABLE gbif AS SELECT * FROM '{gbif_occurence_db_file}' LIMIT 10")




    #Gbif occ to geometry 
    con.execute(f"ALTER TABLE gbif ADD COLUMN geom GEOMETRY")
    con.execute(f"UPDATE gbif SET geom = ST_Point(decimalLongitude, decimalLatitude)")

    #preview_gbif_data(con, "parks")
    #preview_gbif_data(con, "gbif")

    # Create new col for park id 
    #Spatial Join 
    con.execute(f"""CREATE TABLE gbif_joined AS
                    SELECT g.*, p.OBJECTID AS OBJECTID
                    FROM gbif g
                    JOIN parks p
                    ON ST_Within(g.geom, p.geom);""")


    #preview_gbif_data(con, "gbif_joined")
    df = con.execute("SELECT * FROM gbif_joined").df()
    print
    

    print('d')

    

    pass



convert_gbif_csv(gbif_occurence_raw_file, gbif_occurence_db_file)

spatial_join()