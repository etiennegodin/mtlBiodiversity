import duckdb
from pathlib import Path 


db_file_path = Path("data/interim/gbif/gbif_with_parks.parquet")

limit = False
# Create a connection (in-memory or persistent)
con = duckdb.connect() 
 # Install spatial extension 
con.execute("INSTALL spatial;")
con.execute("LOAD spatial;")

if limit:

    con.execute(f"""\
                CREATE TABLE data AS SELECT * FROM '{db_file_path}' 
                WHERE OBJECTID IS NOT NULL
                LIMIT {limit};
                """)
else:
    con.execute(f"""
                CREATE TABLE data AS SELECT * FROM '{db_file_path}' 
                WHERE OBJECTID IS NOT NULL;
                """)


# Species richness per park 
con.execute("""
            CREATE TABLE species_richness AS
            SELECT Nom,
            COUNT(DISTINCT species) AS species_richness
            FROM data
            GROUP BY Nom
            ORDER BY species_richness DESC;
            """)

# Number of observation per year
con.execute("""
            CREATE TABLE annual_observations AS
            SELECT year,
            COUNT(DISTINCT gbifID) AS annual_observations
            FROM data
            GROUP BY year
            ORDER BY year DESC;    
            """)

# Most observed species 
con.execute("""
            
        CREATE TABLE species_observation AS
        SELECT species,
        COUNT(*) AS observations
        FROM data
        GROUP BY species
        ORDER BY observations DESC;
            """)

df = con.execute("SELECT * FROM species_observation").df()
print(df)
# Number of species 
x = con.execute("""
        SELECT COUNT(DISTINCT species) AS total_species
        FROM data;
""").fetchall()
print(x)