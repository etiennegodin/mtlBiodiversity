from pathlib import Path
configfile: "config.yaml"
from math import sqrt

data_dir = Path(config["data_dir"])
raw_dir = data_dir / "raw"
db_dir = data_dir / "db"
scripts_dir = Path("scripts")

raw_shp_path = raw_dir / "geospatial"
int_shp_path = interim_dir / "geospatial"

raw_gbif_path = raw_dir / "gbif" 
int_gbif_path = interim_dir / "gbif" 

db_name = config["duckdb_file"]

# Calculate coordinate uncertainty based on grid size
coord_uncertainty = ((config['grid_size'] * sqrt(2)) / 2)

rule gbif_to_parquet:
    input:
        raw_gbif_path / f"{gbif_raw_file}.csv"
    output:
        int_gbif_path / "gbif_data.parquet"
    params:
        limit = config.get("limit", None)
    script:
        scripts_dir / "01_convertGbifData.py"

rule clean_shapefiles:
    input:
        raw_shp_path/"{name}.shp"
    output:
        int_shp_path/"{name}.shp"
    params:
        crs = config["target_crs"]
    script:
        scripts_dir / "02_clean_shapefiles.py"

rule load_gbif_to_db:
    input:
        int_gbif_path / "gbif_data.parquet"
    output:
        marker = db_dir / ".gbif_table"
    params:
        db_name = config["duckdb_file"],
        marker_file = output.marker,
    script:
        scripts_dir / "03_loadGbifData.py"

rule load_shp_to_db:
    input:
        int_shp_path/"{name}.shp"
    output:
        db_dir/".{name}_table"        
    params:
        db_name = config["duckdb_file"],

    script:
        scripts_dir / "04_load_shapefiles.py"

