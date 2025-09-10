import argparse
from .dataprep.gbif import prep_gbif
from .dataprep.mtl import prep_mtl
from .aggregation.metrics import process_metrics
#from .dashboard.app import run_dashboard

def run_prep(force = False):
    print("Running data prep...")

    prep_mtl(force = force)
    prep_gbif(force = force)

    # Load raw data (CSV, SHP, etc.)
    # Run prep_gbif/prep_mtl
    # Save outputs (e.g. parquet, duckdb, feather)

def run_metrics(force = False):
    print("Running Aggregate Metrics...")
    process_metrics(force = force)

def run_dash():
    print("Launching dashboard...")
    # Load prepped data
    #run_dashboard()

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("command", choices=["prep", "metrics", "dashboard"])
    # optional flag: --force
    parser.add_argument(
        "--force",
        action="store_true",
        help="Force re-run of data prep even if processed files exist"
    )
    args = parser.parse_args()


        
    if args.command == "prep":
        run_prep(force = args.force)
    elif args.command == "metrics":
        run_metrics(force = args.force)
    elif args.command == "dashboard":
        run_dash()
