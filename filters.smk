from pathlib import Path
configfile: "config.yaml"

data_dir = Path(config["data_dir"])
db_dir = data_dir / "db"
scripts_dir = Path("scripts")



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
        wkt_filters = config.get("wkt_filters", [])

    script:
        scripts_dir / "05_filterGbifData.py"