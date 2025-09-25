from pathlib import Path
configfile: "config.yaml"

data_dir = Path(config["data_dir"])
raw_dir = data_dir / "raw"
interim_dir = data_dir / "interim"
processed_dir = data_dir / "processed"
db_dir = data_dir / "db"
scripts_dir = Path("scripts")

raw_shp_path = raw_dir / "geospatial"
int_shp_path = interim_dir / "geospatial"

raw_gbif_path = raw_dir / "gbif" 
int_gbif_path = interim_dir / "gbif" 


# Extract shapefile names from config
shapefile_names = [f.stem for f in raw_shp_path.glob("*.shp")]
gbif_raw_file = [f.stem for f in raw_gbif_path.glob("*.csv")][0]


include: "preprocess.smk"
include: "filters.smk"
include: "analysis.smk"


rule all:
    input:
        expand(db_dir/".{name}_table", name=shapefile_names),
        db_dir / ".gbif_table",
        expand(db_dir/".{name}_sjoin", name=shapefile_names)

