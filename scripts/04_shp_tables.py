from mtlBio.dataprep.duckdb import create_table_from_shp
from mtlBio.core import DuckDBConnection

# Create db 
db = DuckDBConnection(file = snakemake.params.db_name)

for file in snakemake.input:
    create_table_from_shp(file_path= file)
    
    