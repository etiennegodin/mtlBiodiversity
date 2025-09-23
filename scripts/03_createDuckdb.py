from mtlBio.dataprep.spatial_join import create_table_from_shp, create_gbif_table
from mtlBio.core import DuckDBConnection

print('#'*25, '03_createDuckdb', '#'*25,'\n' )


DuckDBConnection.get_connection(file =snakemake.params.db_name)

create_table_from_shp(file_path= snakemake.input[0], limit = snakemake.params.limit)

create_gbif_table(gbif_data_path= snakemake.input[1])


