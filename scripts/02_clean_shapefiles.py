import geopandas as gpd
from mtlBio.dataprep.geospatial import process_shp
print('#'*25, '02_clean_shapefiles', '#'*25,'\n' )

process_shp(input_file = snakemake.input[0], output_file = snakemake.output[0], crs = snakemake.params.crs)
