from mtlBio.preprocess.duckdb import filter_gbif_data
from mtlBio.core import DuckDBConnection

print('#'*25, '05_filterGbifData', '#'*25,'\n' )


# Create db 
db = DuckDBConnection(file = snakemake.params.db_name)


filter_gbif_data(gbif_cols = snakemake.params.gbif_cols,
                    limit = snakemake.params.limit,
                    marker_file = snakemake.params.marker_file,
                    coordUncerFilter = None,
                    taxa_filter = snakemake.params.taxa_filter )
                    
db.conn.close()

