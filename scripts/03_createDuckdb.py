from mtlBio.dataprep.duckdb import create_table_from_shp, create_gbif_table
from mtlBio.core import DuckDBConnection

print('#'*25, '03_createDuckdb', '#'*25,'\n' )


# Create db 
db = DuckDBConnection(file = snakemake.params.db_name)
# Release
create_table_from_shp(file_path= snakemake.input[0])

create_gbif_table(gbif_data_path= snakemake.input[1], limit = snakemake.params.limit )

db.conn.close()

