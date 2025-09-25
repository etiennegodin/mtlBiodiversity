from mtlBio.preprocess.duckdb import create_gbif_table
from mtlBio.core import DuckDBConnection

print('#'*25, '03_loadGbifData', '#'*25,'\n' )


# Create db 
db = DuckDBConnection(file = snakemake.params.db_name)



create_gbif_table(gbif_data_path= snakemake.input[0],
                  marker_file = snakemake.params.marker_file)

db.conn.close()

