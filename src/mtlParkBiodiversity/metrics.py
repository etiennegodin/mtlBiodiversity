import duckdb
from pathlib import Path 
from matplotlib import pyplot as plt
from shapely import wkb,wkt
import geopandas as gpd
from .core import df_inspect

DB_PATH = Path("data/interim/gbif/")
OUTPUT_PATH = Path("data/processed")

# Create a connection (in-memory or persistent)
con = duckdb.connect() 
 # Install spatial extension 
con.execute("INSTALL spatial;")
con.execute("LOAD spatial;")

def save_table(name, geographic_data = False, debug = False):
    print(f"Processing {name} metric")

    df = con.execute(f""""
                     
                     SELECT *,
                     FROM {name}
                     
                     """).df()
    
    df.to_parquet(OUTPUT_PATH / f"{name}.parquet")
    print(f'Saved {name} to parquet file')

    if debug:
        df_inspect(df)

    if geographic_data:
        print("Geographic data, eporting to GeoJSON")
        #Convert park_geom from wkb to shapely geometry
        df["geometry"] = df["park_geom"].apply(lambda x: wkb.loads(bytes(x)))
        gdf = gpd.GeoDataFrame(df, geometry = 'geometry', crs = 'EPSG:4326')
        gdf = gdf.drop(columns = ['park_geom'])
        
        #Clean Nan values 
        gdf = gdf.dropna(how="any").reset_index(drop=True)
        # Inspect
        if debug:
            df_inspect(gdf)

        gdf.to_file(OUTPUT_PATH / f"{name}.geojson", driver = 'GeoJSON')

# Species richness per park 
park_metrics_query = """
CREATE TABLE parks AS
SELECT park_name,
COUNT(DISTINCT gbifID) AS observation, -- dsda-- 
COUNT(DISTINCT kingdom) AS kingdom_richness,
COUNT(DISTINCT family) AS family_richness,
COUNT(DISTINCT genus) AS genus_richness,
COUNT(DISTINCT species) AS species_richness, 

FROM data
WHERE park_name IS NOT NULL
GROUP BY park_name
ORDER BY species_richness DESC

LIMIT 1000;

"""
shannon_index_query = """

CREATE OR REPLACE VIEW park_shannon_index AS
WITH species_count AS
(SELECT park_name, species, 
COUNT(*) AS species_count,
FROM data
WHERE park_name IS NOT NULL
GROUP BY park_name, species,
ORDER BY park_name
),

park_totals AS(

SELECT park_name,
SUM(species_count) AS total_count
FROM species_count
GROUP BY park_name,
),

proportions AS(

SELECT s.park_name,
        s.species,
        s.species_count,
        t.total_count,
        s.species_count * 1.0 / ( t.total_count *1.0 )AS p
        
FROM species_count s
JOIN park_totals t
ON s.park_name = t.park_name
)


SELECT park_name,

    -SUM(p*LN(p)) AS shannon_index
FROM proportions
GROUP BY park_name;

"""

join_shannon = """
SELECT p.*, s.shannon_index
FROM parks p
LEFT JOIN park_shannon_index s
on p.park_name = s.park_name

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

    if test:
        db_file_path = DB_PATH / "_test_gbif_with_parks.parquet"
    else:
        db_file_path = DB_PATH / "gbif_with_parks.parquet"



    print("Loading gbif with parks data")

    if test:
        print(f"Loading data with limit set to {limit}")
        con.execute(f"""
                    CREATE TABLE data AS SELECT * FROM '{db_file_path}' 
                    WHERE park_id IS NOT NULL
                    LIMIT {limit};
                    """)
    else:
        con.execute(f"""
                        CREATE TABLE data AS SELECT * FROM '{db_file_path}' 
                        WHERE park_id IS NOT NULL;
                        """)
        

    con.execute(park_metrics_query)
    df = con.execute("""SELECT *, FROM parks""").df()
    con.execute(shannon_index_query)

    df = con.execute(join_shannon).df()
    print(df)

        

    #save_table('parks', geographic_data= True , debug=False)
    #save_table('annual_observations', query = annual_observations_query)
    #save_table('most_observed_species', query = most_observed_species_query)
    #save_table('species_count', query = species_count_query)

    con.close()
