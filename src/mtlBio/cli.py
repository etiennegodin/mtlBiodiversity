import argparse
import subprocess
from pathlib import Path
import sys
from mtlBio.core  import raw_data_read_only
from mtlBio.dataprep.gbif import prep_gbif_data, gbif_spatial_joins
from mtlBio.dataprep.geospatial import prep_geospatial
from mtlBio.metrics.main import process_all_metrics



def run_prep(force = False, test = False, limit = None, colab = False):
    print("Running data prep...")

    if force:
        user_force = input("You have chosen to force re-run GBIF data prep. This may take a while. Do you want to continue? (y/n): ")
        if user_force.lower() == 'y':
            force = True
        else:
            force = False
            print("Skipping GBIF data prep.")

    # Load raw data (CSV, SHP, etc.)
    # Run prep_gbif/prep_mtl
    # Save outputs (e.g. parquet, duckdb, feather)
def run_geospatial(force = False):
    prep_geospatial(force = force)

def run_gbif(force = False, test = False, limit = None, colab = False):
    
    gbif_parquet = prep_gbif_data(force =force, test =test, limit =limit)

    gbif_spatial_joins(gbif_occurence_db_file = gbif_parquet, force =force, test =test, limit =limit)

def run_metrics(force = False, test = False, limit = None):
    print("Running Aggregate Metrics...")
    process_all_metrics(force = force, test = test, limit = limit)

def run_app():
    print("Launching dashboard...")
    # Load prepped data
    #run_appboard()
    app_folder = Path(__file__).parent.parent.parent
    # Run Streamlit app with streamlit subprocess
    subprocess.run([sys.executable, "-m", "streamlit", "run", f"{app_folder}/streamlit_app.py"])
    
    
def launch_debugger():
    
    import debugpy
    debugpy.listen(("localhost", 5678))
    print("Waiting for debugger attach...")
    debugpy.wait_for_client()  # Pause until VS Code connects
    
def main():
    raw_data_read_only('data/raw', debug = False)
    parser = argparse.ArgumentParser()
    parser.add_argument("command", choices=["geo", "gbif", 'sjoin', "metrics", "app", 'full'])
    # optional flag: --force
    parser.add_argument(
        "--force",
        action="store_true",
        help="Force re-run of data prep even if processed files exist"
    )
    parser.add_argument("--debug", action= 'store_true', help = 'Run debugger')
    parser.add_argument("--test", action= 'store_true', help = 'Run in test mode')
    parser.add_argument("--limit", type= int, help = 'Limit the number of rows processed (for testing purposes)')

    args = parser.parse_args()
    if args.debug:
        launch_debugger()
    # Default limit to 1000 if test mode is on and limit is not specified
    if args.limit is None:
        args.limit = 1000
    if args.command == 'full':
        run_geospatial(force = args.force)
        run_gbif(force = args.force, test = args.test, limit = args.limit)
        run_metrics(force = args.force, test = args.test, limit = args.limit)
        run_app()
    elif args.command == "geo":
        run_geospatial(force = args.force)
    elif args.command == "gbif":
        run_gbif(force = args.force, test = args.test, limit = args.limit)
    elif args.command == "metrics":
        run_metrics(force = args.force, test = args.test, limit = args.limit)
    elif args.command == "app":
        run_app()
        
        

if __name__ == "__main__":
    main()