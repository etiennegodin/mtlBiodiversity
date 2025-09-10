import argparse
from .dataprep.gbif import prep_gbif
from .dataprep.mtl import prep_mtl
from .metrics import process_metrics
#from .dashboard.app import run_dashboard

def run_prep(force = False, test = False, limit = None):
    print("Running data prep...")

    prep_mtl(force = force, )
    prep_gbif(force = force, test = test, limit = limit)


    # Load raw data (CSV, SHP, etc.)
    # Run prep_gbif/prep_mtl
    # Save outputs (e.g. parquet, duckdb, feather)

def run_metrics(force = False, test = False, limit = None):
    print("Running Aggregate Metrics...")
    process_metrics(force = force, test = test, limit = limit)

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
    parser.add_argument("--test", action= 'store_true', help = 'Run in test mode')
    parser.add_argument("--limit", type= int, help = 'Limit the number of rows processed (for testing purposes)')
    args = parser.parse_args()

    # Default limit to 10000 if test mode is on and limit is not specified
    if args.limit is None:
        args.limit = 10000
        
    if args.command == "prep":
        run_prep(force = args.force, test = args.test, limit = args.limit)
    elif args.command == "metrics":
        run_metrics(force = args.force, test = args.test, limit = args.limit)
    elif args.command == "dashboard":
        run_dash()
