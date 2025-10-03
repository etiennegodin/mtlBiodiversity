from pathlib import Path
configfile: "config.yaml"

data_dir = Path(config["data_dir"])
processed_dir = data_dir / "processed"
db_dir = data_dir / "db"
r_scripts = Path("R")

rule RsamplingBias:
    input:
        db_dir/".grid_sjoin",
        "renv.lock"
    output:
        r_scripts/"ids_samplingBias.rds"
    resources:
        db = 1       
    shell:
        """
        Rscript R/samplingBias.R
        """

rule Ranalysis:
    input:
        r_scripts/"ids_samplingBias.rds"
    output:
        out = processed_dir/".{name}.parquet"
    resources:
        db = 1     
    params:
        db_name = config["duckdb_file"],
        table = "{name}_sjoin",
        group_col = "{name}_id"        
    conda:
        "envs/r_env.yaml"
    shell:
        "Rscript R/group_analysis.R"


