#Snakefile
configfile: "config.yaml"

SAMPLES = list(config["geo_samples"].keys())

rule all:
    input:
        expand("data/interim/geospatial/{sample}_inter.shp", sample=SAMPLES)


rule clean_geospatial:
    input:
        lambda wildcards: config["geo_samples"][wildcards.sample]["file"]
    output:
        "data/interim/geospatial/{sample}_inter.shp"
    script:
        "scripts/01_clean_geospatial.py"