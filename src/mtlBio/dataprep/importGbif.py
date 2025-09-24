from mtlBio.core import convertToPath
import duckdb

def main(input_path, output_path, limit = None):

    if limit == 'None':
        limit = None
        
    query = None
    
    input_path = convertToPath(input_path)
    output_path = convertToPath(output_path)
    print(f'Converting to {input_path} to {output_path} file with limit set to {limit} ')
    query = f"COPY (SELECT * FROM '{input_path}' {'LIMIT ' + str(limit) if (limit is not None) else ''}) TO '{output_path}' (FORMAT PARQUET)"
    if query is not None:
        try:
            duckdb.query(query)
            print('-'*25,'Successfuly convert gbif csv data to parquet','-'*25, '\n')
        except Exception as e:
            print(f"Failed to convert csv to parquet: {e}")
            return False
