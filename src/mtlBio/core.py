import geopandas as gpd
from pathlib import Path
import pandas as pd
from shapely.geometry import Polygon, MultiPolygon, Point
import pathlib
import stat
import duckdb
import tkinter as tk 
from tkinter import filedialog
from jinja2 import Template

SQL_PATH = Path("data/sql_queries")

class DuckDBConnection:
    _instance = None

    @classmethod
    def get_connection(cls, file):
        if cls._instance is None:
            cls._instance = duckdb.connect(file)
        return cls._instance


def convertToPath(file_path: str):
    path = Path(file_path)
    return path
    
def find_files(folder_path: Path, expected_files : list = None, suffix = None, debug :bool= False):
    files_found = []
    files_not_found = []

    
    #Check if suffix is provided
    if suffix is not None:
        # If suffix is only a string format to list
        if isinstance(suffix, str):
            suffix_list = [suffix]
        else:
            suffix_list = suffix
            
        # Check if written with a dot - fix if not
        final_suffix  = []
        for s in suffix_list:
            if not s.startswith('.'):
                s = f'.{s}'
                final_suffix.append(s)
            else:
                final_suffix.append(s)
        print(final_suffix)
        files = [f for f in folder_path.rglob("*") if f.suffix in final_suffix]
        
        print(f"Found {len(files)} file in {folder_path}")
        for e_file in expected_files:
            e_file_found = False
            for f in files:
                f_name = f.stem.lower()
                if debug:
                    print(f_name)
                    print(e_file)
                if e_file in f_name:
                    files_found.append(f)
                    e_file_found = True
                    break
                
            if not e_file_found: 
                files_not_found.append(e_file)
                    
        print(files_found)
        for file in files_not_found:
            print(f'{file} not found in {folder_path}, please provide')
            file = select_file()   
            files_found.append(file)
    
    if len(files_found) == len(expected_files):
        return files_found



def read_sql_template(file_name : str = None):
    global SQL_PATH
    with open(SQL_PATH / f'{file_name}.sql', 'r') as f:
        sql_template = Template(f.read())
        
    return sql_template

def read_sql_file(file_name):
    global SQL_PATH
    with open(SQL_PATH / f'{file_name}.sql', 'r') as f:
        sql_query = f.read()
    
    return sql_query

def select_file():
    root = tk.Tk()
    root.withdraw()  # hide the root window
    file_path = filedialog.askopenfilename(
        title="Select a file",
        filetypes=[("All files", "*.*"), ("Text files", "*.txt"), ("CSV files", "*.csv")]
    )
    return file_path

def convert_df_to_gdf(df : pd.DataFrame, lat_col : str = 'decimalLatitude', long_col : str = 'decimalLongitude', crs = 4326, verbose = False):
    gdf = gpd.GeoDataFrame(df, geometry=[Point(xy) for xy in zip(df["decimalLongitude"], df["decimalLatitude"])] , crs = 4326 )
    return gdf

def raw_data_read_only(path, debug = True):
    # Set read-only permissions for all files in the specified directory
    
    folder = pathlib.Path(path)

    for file in folder.rglob("*"):
        if file.is_file():
            try:    
                file.chmod(stat.S_IREAD)
            except Exception as e:
                if debug:
                    print(f"Failed to change permissions for {file}: {e}")

def df_inspect(df : pd.DataFrame, limit = 100):
    print('-'*25, 'Preview', '-'*25)
    print(df.head(limit))
    print('/'*25, 'Columns', '/'*25)
    print(df.columns)
    print('/'*25, 'Dtypes', '/'*25)
    print(df.dtypes)
    print('/'*25, 'Shape', '/'*25)
    print(df.shape)
    print('-'*25, 'Null values', '-'*25)
    print(df.isnull().sum())
    print('/'*50)

def convert_crs(gdf : gpd.GeoDataFrame, target_crs = 4326, verbose = False):
    """
    Converts gdf to taregt crs 
    """
    crs_init = gdf.crs

    if crs_init != target_crs:
        try:
            gdf = gdf.to_crs(target_crs)
            crs_new = gdf.crs
            if crs_init != crs_new:
                return gdf
            elif crs_init == crs_new:
                print('Failed to convert crs')
                return None

        except Exception as e:
            print(f"Failed to convert gdf to taregt crs {e}")
            return None

    elif crs_init == target_crs:
        if verbose:
            print(f"Gdf already with target crs {target_crs}")
        return gdf



def check_common_crs(gdf1 : gpd.GeoDataFrame, gdf2 : gpd.GeoDataFrame):
    """
    Check if 2 geodataframes have the same crs
    """
    if gdf1.crs == gdf2.crs:
        return True
    else:
        return False

def clip_to_region(input : gpd.GeoDataFrame, clip_region : gpd.GeoDataFrame):
    """
    Clip input gdf by other

    Checks crs compatibility and auto convert crs if not the same 
    """
    if not check_common_crs(input, clip_region):
        target_crs = input.crs
        clip_region = convert_crs(clip_region, target_crs)
    try:
        clip_region_merged = clip_region.union_all()
        clipped_input = input.clip(clip_region_merged)
        return clipped_input 

    except Exception as e:
        print(f"Failed to clip input by clip_region : {e}")

    
