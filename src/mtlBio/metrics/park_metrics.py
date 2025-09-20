import duckdb
from pathlib import Path 
from matplotlib import pyplot as plt
from shapely import wkb
import geopandas as gpd
from mtlBio.core import df_inspect

DB_PATH = Path("data/interim/gbif/")
OUTPUT_PATH = Path("data/processed")
SQL_PATH = Path("data/sql_queries")

def read_sql_file(file_name):
    with open(SQL_PATH / f'{file_name}.sql', 'r') as f:
        park_metrics_sql = f.read()
    
    return park_metrics_sql


def save_table(name, geographic_data = False, debug = False, con = None):
    print(f"Processing {name} metric")

    
# Save table to Parquet
    con.execute(f"COPY {name} TO '{OUTPUT_PATH}/{name}.parquet' (FORMAT PARQUET);")
    print(f'Saved {name} metric to parquet file')

    if geographic_data:
        print("Geographic data, exporting to GeoJSON")

        df = con.execute(f"SELECT * FROM {name}").df()
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

CREATE OR REPLACE TABLE parks AS

WITH kingdom_richness AS(
SELECT park_name, kingdom,
COUNT(*) AS count
FROM data
GROUP BY park_name, kingdom,

),

class_richness AS(
SELECT park_name, class,
COUNT(*) AS count
FROM data
GROUP BY park_name, class,
),

taxa_richness AS(
SELECT k.count AS kingdom_count,
SELECT c.count AS class_count,

FROM kingdom_richness k
JOIN class_richness c
ON k.park_name = c.park_name

)


SELECT park_name,
COUNT(DISTINCT gbifID) AS observation, -- dsda-- 
COUNT(DISTINCT kingdom) AS kingdom_richness,
COUNT(DISTINCT family) AS family_richness,
COUNT(DISTINCT genus) AS genus_richness,
COUNT(DISTINCT species) AS species_richness,

ANY_VALUE(ST_AsWKB(park_geom)) AS park_geom,

t.kingdom_count,
t.class_count

FROM data
JOIN taxa_richness t
ON data.park_name = t.park_name
WHERE park_name IS NOT NULL
GROUP BY park_name
ORDER BY species_richness DESC

;

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


def park_metrics(force = False, test = False, limit = None):

    if test:
        db_file_path = DB_PATH / "_test_gbif_with_parks.parquet"
    else:
        db_file_path = DB_PATH / "gbif_with_parks.parquet"

    # Create a connection (in-memory or persistent)
    con = duckdb.connect('data/db/parks_metrics.duckdb') 



    print("Loading gbif with parks data")

    if test:
        print(f"Loading data with limit set to {limit}")
        con.execute(f"""
                    CREATE OR REPLACE TABLE data AS
                    SELECT * FROM '{db_file_path}' 
                    WHERE park_id IS NOT NULL
                    LIMIT {limit};
                    """)
    else:
        con.execute(f"""
                    CREATE OR REPLACE TABLE data AS
                    SELECT * FROM '{db_file_path}' 
                    WHERE park_id IS NOT NULL
                    ;
                    """)
        

    # Install spatial extension 
    con.execute("INSTALL spatial;")
    con.execute("LOAD spatial;")

    park_metrics_sql = read_sql_file('parks')
    con.execute(park_metrics_sql)

 

    shannon_index_sql = read_sql_file('shannon_index')
    con.execute(shannon_index_sql)



    con.execute("ALTER TABLE parks ADD COLUMN shannon_index DOUBLE")
    con.execute("""
        UPDATE parks
        SET shannon_index = s.shannon_index
        FROM park_shannon_index s
        WHERE parks.park_name = s.park_name;
    """)

    df = con.execute("SELECT * FROM parks").df()
    print(df)

    save_table('parks', geographic_data= True , debug=False, con = con)

    #taxa_group_sql = read_sql_file('park_taxa')
    #con.execute(taxa_group_sql)
    #save_table('park_taxa', geographic_data= True , debug=False, con = con)


    con.close()
