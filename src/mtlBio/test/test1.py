

import duckdb
from shapely import wkb, wkt
import geopandas as gpd

file = 'data\interim\joined\occurences_quartiers.parquet'
con = duckdb.connect()

df = con.execute(f"""SELECT * FROM '{file}' WHERE geom is NOT NULL""").df()
print(df.head())
#Convert park_geom from wkb to shapely geometry
df["geometry"] = df["geom"].apply(lambda x: wkb.loads(bytes(x)))
gdf = gpd.GeoDataFrame(df, geometry = 'geometry', crs = 'EPSG:4326')

gdf.to_file('src/mtlBio/test/test.shp')    