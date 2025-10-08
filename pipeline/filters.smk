from pathlib import Path
configfile: "config.yaml"

data_dir = Path(config["data_dir"])
db_dir = Path(config["db_dir"])
scripts_dir = Path(config["scripts_dir"])

rule filter_gbif_data:
    input:
        db_dir / ".gbif_table"
    output:
        marker = db_dir / ".filtered_gbif_table"
    params:
        gbif_cols = config["gbif_columns"],
        limit = config.get("limit", None),
        db_name = config["duckdb_file"],
        marker_file = output.marker,
        coordUncerFilter = coord_uncertainty,
        taxa_filter =  config["taxa_filter"]

    script:
        scripts_dir / "05_filterGbifData.py"

