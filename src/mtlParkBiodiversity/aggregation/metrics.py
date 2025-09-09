import duckdb
from pathlib import Path 
from matplotlib import pyplot as plt

db_file_path = Path("data/interim/gbif/gbif_with_parks.parquet")
output_path = Path("data/processed")
limit = False

# Create a connection (in-memory or persistent)
con = duckdb.connect() 
 # Install spatial extension 
con.execute("INSTALL spatial;")
con.execute("LOAD spatial;")

def create_metric(name, query = None, save = False, limit = None):
    df = con.execute(query).df()

    if save:
        df.to_parquet(output_path / f"{name}.parquet")


if limit:

    con.execute(f"""\
                CREATE TABLE data AS SELECT * FROM '{db_file_path}' 
                WHERE OBJECTID IS NOT NULL
                LIMIT 10000;
                """)
else:
    con.execute(f"""
                CREATE TABLE data AS SELECT * FROM '{db_file_path}' 
                WHERE OBJECTID IS NOT NULL;
                """)



# Species richness per park 
species_richness_query = """
                            SELECT Nom,
                            COUNT(DISTINCT species) AS species_richness
                            FROM data
                            GROUP BY Nom
                            ORDER BY species_richness DESC;
                            """

create_metric('species_richness', query= species_richness_query, save= True )


# Number of observation per year
annual_observations_query = """
                            SELECT year,
                            COUNT(DISTINCT gbifID) AS annual_observations
                            FROM data
                            GROUP BY year
                            ORDER BY year DESC;    
                            """
create_metric('annual_observations', query = annual_observations_query, save = True)

# Most observed species 

most_observed_species_query = """            
                            SELECT species,
                            COUNT(*) AS observations
                            FROM data
                            GROUP BY species
                            ORDER BY observations DESC;
                            """

create_metric('most_observed_species', query = most_observed_species_query, save = True)

# Number of species 
species_count_query = """
        SELECT COUNT(DISTINCT species) AS total_species
        FROM data;
"""

create_metric('species_count', query = species_count_query, save = True)