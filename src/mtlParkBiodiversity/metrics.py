import duckdb
from pathlib import Path 
from matplotlib import pyplot as plt
from shapely import wkb,wkt
import geopandas as gpd
from .core import df_inspect

db_file_path = Path("data/interim/gbif/gbif_with_parks.parquet")
output_path = Path("data/processed")
limit = False

# Create a connection (in-memory or persistent)
con = duckdb.connect() 
 # Install spatial extension 
con.execute("INSTALL spatial;")
con.execute("LOAD spatial;")

def create_metric(name, query = None, file_type = "parquet", debug = False):
    print(f"Processing {name} metric")
    df = con.execute(query).df()

    if file_type == "parquet":

        # Inspect
        if debug:
            df_inspect(df)
        df.to_parquet(output_path / f"{name}.{file_type}")

    elif file_type == 'geojson':
        #Convert park_geom from wkb to shapely geometry
        print("Converting park_geom from wkb to shapely geometry")
        df["geometry"] = df["park_geom"].apply(lambda x: wkb.loads(bytes(x)))
        gdf = gpd.GeoDataFrame(df, geometry = 'geometry', crs = 'EPSG:4326')
        gdf = gdf.drop(columns = ['park_geom'])
        
        #Clean Nan values 
        gdf = gdf.dropna(how="any").reset_index(drop=True)
        # Inspect
        if debug:
            df_inspect(gdf)

        gdf.to_file(output_path / f"{name}.{file_type}", driver = 'GeoJSON')

# Species richness per park 
species_richness_query = """
                            SELECT Nom,
                            ST_AsWKB(park_geom) AS park_geom,
                            COUNT(DISTINCT species) AS species_richness
                            FROM data
                            GROUP BY Nom, park_geom
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

#Shannon index
# Type species per park 
# genus count per park 
# family count per park
# order count per park



def process_metrics(force = False, test = False, limit = None):

    print("Loading gbif with parks data")

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
        

    create_metric('species_richness', query= species_richness_query, file_type = 'geojson', debug=True)
    #create_metric('annual_observations', query = annual_observations_query)
    #create_metric('most_observed_species', query = most_observed_species_query)
    #create_metric('species_count', query = species_count_query)

    con.close()
