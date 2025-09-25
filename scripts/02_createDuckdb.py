from mtlBio.preprocess.duckdb import create_gbif_table
from mtlBio.core import DuckDBConnection

print('#'*25, '03_createDuckdb', '#'*25,'\n' )


# Create db 
db = DuckDBConnection(file = snakemake.params.db_name)



create_gbif_table(gbif_data_path= snakemake.input[0],
                  limit = snakemake.params.limit,
                  marker_file = snakemake.params.marker_file,
                  coordUncerFilter = snakemake.params.coordUncerFilter )

db.conn.close()

