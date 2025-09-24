from pathlib import Path
configfile: "config.yaml"

data_dir = Path(config["data_dir"])
raw_dir = data_dir / "raw"
interim_dir = data_dir / "interim"
processed_dir = data_dir / "processed"
db_dir = data_dir / "db"

raw_shp_path = raw_dir / "geospatial"
int_shp_path = interim_dir / "geospatial"

raw_gbif_path = raw_dir / "gbif" 
int_gbif_path = interim_dir / "gbif" 

scripts_dir = Path("scripts")

# Extract shapefile names from config
shapefile_names = [f.stem for f in raw_shp_path.glob("*.shp")]
shapefile_paths = [f for f in raw_shp_path.glob("*.shp")]

gbif_raw_file = [f.stem for f in raw_gbif_path.glob("*.csv")][0]


rule all:
    input:
        expand(db_dir/".{name}_done", name=shapefile_names),
        db_dir / ".gbif_done"

rule load_gbif_data:
    input:
        raw_gbif_path / f"{gbif_raw_file}.csv"
    output:
        int_gbif_path / "gbif_data.parquet"
    params:
        limit = config.get("limit", None)
    script:
        scripts_dir / "01_loadGbifData.py"

rule clean_shapefiles:
    input:
        raw_shp_path/"{name}.shp"
    output:
        int_shp_path/"{name}.shp"
    params:
        crs = config["target_crs"]
    script:
        scripts_dir / "02_clean_shapefiles.py"

rule create_duckdb:
    input:
        int_gbif_path / "gbif_data.parquet"
    output:
        marker = db_dir / ".gbif_done"
    params:
        limit = config.get("limit", None),
        db_name = config["duckdb_file"],
        marker_file = output.marker

    script:
        scripts_dir / "03_createDuckdb.py"

rule create_shp_tables:
    input:
        int_shp_path/"{name}.shp"
    output:
        db_dir/".{name}_done"
    params:
        db_name = config["duckdb_file"],

    script:
        scripts_dir / "04_shp_tables.py"

rule test:
    input:
        db_dir/".{name}_done"
    output:
        config["duckdb_file"]
# make .ready files of tables instaed of updating main .duckdb as snakefile doesnt understand and skips step 4 
"""
rule spatial_joins:
    input:
        expand(int_shp_path/"{shapefile}", shapefile=shapefile_names)
    output:
        config["duckdb_file"]
    script:
        scripts_dir / "04_spatial_join.py"

"""