from pathlib import Path
configfile: "config.yaml"

data_dir = Path(config["data_dir"])
raw_dir = Path(config["raw_dir"])
interim_dir : Path(config["interim_dir"])
db_dir = Path(config["db_dir"])
scripts_dir = Path(config["scripts_dir"])

raw_shp_path = Path(config["raw_shp_path"])
int_shp_path = Path(config["int_shp_path"])
raw_gbif_path = Path(config["raw_gbif_path"])
int_gbif_path = Path(config["int_gbif_path"])

db_name = Path(config["duckdb_file"])


# Extract shapefile names from config
shapefile_names = [f.stem for f in raw_shp_path.glob("*.shp")]
gbif_raw_file = [f.stem for f in raw_gbif_path.glob("*.csv")][0]


include: "pipeline/preprocess.smk"
include: "pipeline/filters.smk"
include: "pipeline/analysis.smk"


rule all:
    input:
        expand(db_dir/".{name}_table", name=shapefile_names),
        db_dir / ".gbif_table",
        expand(db_dir/".{name}_sjoin", name=shapefile_names)

