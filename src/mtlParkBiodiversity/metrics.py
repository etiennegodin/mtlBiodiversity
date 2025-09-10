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

def create_metric(name, query = None):
    print(f"Processing {name} metric")
    df = con.execute(query).df()
    df.to_parquet(output_path / f"{name}.parquet")

# Species richness per park 
species_richness_query = """
                            SELECT Nom,
                            COUNT(DISTINCT species) AS species_richness
                            FROM data
                            GROUP BY Nom
                            ORDER BY species_richness DESC;
                            """

# Number of observation per year
annual_observations_query = """
                            SELECT year,
                            COUNT(DISTINCT gbifID) AS annual_observations
                            FROM data
                            GROUP BY year
                            ORDER BY year DESC;    
                            """

# Most observed species 
most_observed_species_query = """            
                            SELECT species,
                            COUNT(*) AS observations
                            FROM data
                            GROUP BY species
                            ORDER BY observations DESC;
                            """

# Number of species 
species_count_query = """
        SELECT COUNT(DISTINCT species) AS total_species
        FROM data;
"""

def process_metrics(force = False, test = False, limit = None):

    print("Loading gbif with parks data")
    if limit is None:
        limit = 1000

    if test:
        print(f"Loading data with limit set to {limit}")
        con.execute(f"""
                    CREATE TABLE data AS SELECT * FROM '{db_file_path}' 
                    WHERE OBJECTID IS NOT NULL
                    LIMIT {limit};
                    """)
    else:
        con.execute(f"""
                        CREATE TABLE data AS SELECT * FROM '{db_file_path}' 
                        WHERE OBJECTID IS NOT NULL;
                        """)
        

    create_metric('species_richness', query= species_richness_query )
    create_metric('annual_observations', query = annual_observations_query)
    create_metric('most_observed_species', query = most_observed_species_query)
    create_metric('species_count', query = species_count_query)

    con.close()
