#Snakefile
configfile: "config.yaml"

SAMPLES = list(config["geo_samples"].keys())

rule all:
    input:
        expand("data/interim/geospatial/{sample}_inter.shp", sample=SAMPLES),
        "data/interim/gbif/gbif_data.parquet"

rule clean_geospatial:
    input:
        lambda wildcards: config["geo_samples"][wildcards.sample]["file"]
    output:
        "data/interim/geospatial/{sample}_inter.shp"
    script:
        "scripts/01_clean_geospatial.py"

rule load_gbif_data:
    input:
        "data/raw/gbif/gbif_occurences.csv"
    output:
        "data/interim/gbif/gbif_data.parquet"
    params:
        limit = config.get("limit", None)
    script:
        "scripts/02_loadGbifData.py"

