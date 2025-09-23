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

print(gbif_raw_file)
print(shapefile_names)
print(int_shp_path)

rule all:
    input:
        config["duckdb_file"]

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
        shapefile=lambda wildcards: f"{raw_shp_path}/{wildcards.shapefile}.shp"
    output:
        cleaned=int_shp_path/"{shapefile}_clean.shp"
    params:
        crs = config["target_crs"]
    script:
        scripts_dir / "02_clean_shapefiles.py"

rule create_duckdb:
    input:
        expand(int_shp_path/"{shapefile}_clean.shp", shapefile=shapefile_names),
        int_gbif_path / "gbif_data.parquet"
    output:
        config["duckdb_file"]
    params:
        limit = config.get("limit", None),
        db_name = config["duckdb_file"]
    script:
        scripts_dir / "03_createDuckdb.py"

