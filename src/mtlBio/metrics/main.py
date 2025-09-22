from mtlBio.metrics import old_park_metrics
from mtlBio.core import DuckDBConnection, read_sql_template

def process_all_metrics(force = False, test = False, limit = None):

    con = DuckDBConnection.get_connection()
    # Run park metrics

    # Install spatial extension 
    con.execute("INSTALL spatial;")
    con.execute("LOAD spatial;")
    
    taxa_richness_template = read_sql_template('taxa_richness')
    ``
    taxa_richness_query_grid = taxa_richness_template.render(
                    new_table_name = 'grid_taxa_richness',
                    group_by = 'id',
                    table_name = '_test_occurences_grid'            
    )
    
    try:
        con.execute(taxa_richness_query_grid)
    except Exception as e:
        print(e)
    
    """
    
    
    file = gbif_spatial_join
    con = duckdb.connect(x)
    
    run aggregate on tables by specifying groupby column ( park_id, id, qr_id)
    
    save to arquet + gson 
    
    """
    
    
if __name__ == "__main__":
    process_all_metrics()