from pathlib import Path
configfile: "config.yaml"

data_dir = Path(config["data_dir"])
raw_dir = data_dir / "raw"
interim_dir = data_dir / "interim"
processed_dir = data_dir / "processed"
db_dir = data_dir / "db"
scripts_dir = Path("scripts")





rule grid_spatial_join:
    input:
        db_dir / ".filtered_gbif_table",
        db_dir/".grid_table"
    output:
        db_dir/".grid_sjoin"
    resources:
        db_con =1 
    group: "duckdb"
        
    params:
        db_name = config["duckdb_file"],        
    script:
        scripts_dir / "05_grid_sjoin.py"

rule shp_spatial_join:
    input:
        db_dir/".grid_sjoin",
        db_dir/".{name}_table",
        db_lock
    group: "duckdb"

    output:
        db_dir/".{name}_sjoin"
    resources:
        db_con =1         
    params:
        db_name = config["duckdb_file"],        
    script:
        scripts_dir / "06_shp_sjoin.py"
