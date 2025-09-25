from mtlBio.analysis.spatial import grid_spatial_join
from mtlBio.core import DuckDBConnection

print('#'*25, '06_shp_sjoin', '#'*25,'\n' )

try:
    print('-'*100, 'connection')
    db = DuckDBConnection(file = snakemake.params.db_name)
except Exception as e:
    print('*'*100,'fail shp tble:', e)

# Transform input marker file back to table name 
left_table_name = snakemake.input[0].split(sep=".")[-1] # Keeping grid_sjoin table
right_table_name = snakemake.input[1].split(sep=".")[-1].split(sep='_')[0]

#Join grid first
grid_spatial_join(left_table_name = left_table_name, right_table_name= right_table_name, marker_file = snakemake.output[0])


db.conn.close()
 
    