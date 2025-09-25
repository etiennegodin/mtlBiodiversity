import argparse
import subprocess
from pathlib import Path
import sys
from mtlBio.config import configs

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

    if args.command == "app":
        run_app()
        
        

if __name__ == "__main__":
    main()