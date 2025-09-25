from mtlBio.analysis.spatial import grid_spatial_join
from mtlBio.core import DuckDBConnection

print('#'*25, '05_grid_spatial_join', '#'*25,'\n' )


db = DuckDBConnection(file = snakemake.params.db_name)

# Transform input marker file back to table name 


left_table_name = snakemake.input[0].split(sep=".")[-1].split(sep='_')[0]
right_table_name = snakemake.input[1].split(sep=".")[-1].split(sep='_')[0]

print(left_table_name)
print(right_table_name)

#Join grid first
grid_spatial_join(left_table_name = left_table_name, right_table_name= right_table_name, marker_file = snakemake.output[0])


db.conn.close()
 
    