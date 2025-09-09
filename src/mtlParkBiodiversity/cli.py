import argparse
#from .dataprep.gbif import prep_gbif
#from .dataprep.mtl import prep_mtl
#from .dashboard.app import run_dashboard

def run_prep():
    print("Running data prep...")
    # Load raw data (CSV, SHP, etc.)
    # Run prep_gbif/prep_mtl
    # Save outputs (e.g. parquet, duckdb, feather)

def run_dash():
    print("Launching dashboard...")
    # Load prepped data
    #run_dashboard()

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("command", choices=["prep", "dashboard"])
    args = parser.parse_args()

    if args.command == "prep":
        run_prep()
    elif args.command == "dashboard":
        run_dash()
