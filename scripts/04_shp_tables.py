from mtlBio.dataprep.duckdb import create_table_from_shp
from mtlBio.core import DuckDBConnection

print('#'*25, '04_shp_tables', '#'*25,'\n' )

# Create db 
db = DuckDBConnection(file = snakemake.params.db_name)

create_table_from_shp(file_path= snakemake.input[0], marker_file = snakemake.output[0])
    
db.conn.close()
