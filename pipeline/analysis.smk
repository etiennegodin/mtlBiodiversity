from pathlib import Path
configfile: "config.yaml"

data_dir = Path(config["data_dir"])
raw_dir = Path(config["raw_dir"])
db_dir = Path(config["db_dir"])
scripts_dir = Path(config["scripts_dir"])
raw_shp_path = Path(config["raw_shp_path"])
int_shp_path = Path(config["int_shp_path"])
raw_gbif_path = Path(config["raw_gbif_path"])
int_gbif_path = Path(config["int_gbif_path"])
db_name = Path(config["duckdb_file"])


rule grid_spatial_join:
    input:
        db_dir / ".filtered_gbif_table",
        db_dir/".grid_table"
    output:
        db_dir/".grid_sjoin"
    params:
        db_name = config["duckdb_file"]        
    script:
        scripts_dir / "06_grid_sjoin.py"

rule shp_spatial_join:
    input:
        db_dir/".grid_sjoin",
        db_dir/".{name}_table"
    output:
        db_dir/".{name}_sjoin"
    resources:
        db = 1       
    params:
        db_name = config["duckdb_file"]        
    script:
        scripts_dir / "07_shp_sjoin.py"
