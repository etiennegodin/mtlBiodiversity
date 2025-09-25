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
    params:
        db_name = config["duckdb_file"],        
    script:
        scripts_dir / "06_grid_sjoin.py"

rule shp_spatial_join:
    input:
        db_dir/".grid_sjoin",
        db_dir/".{name}_table"
    output:
        db_dir/".{name}_sjoin"       
    params:
        db_name = config["duckdb_file"],        
    script:
        scripts_dir / "07_shp_sjoin.py"
