import geopandas as gpd
import pandas as pd
from shapely.geometry import Polygon, MultiPolygon, Point

def convert_df_to_gdf(df : pd.DataFrame, lat_col : str = 'decimalLatitude', long_col : str = 'decimalLongitude', crs = 4326, verbose = False):
    gdf = gpd.GeoDataFrame(df, geometry=[Point(xy) for xy in zip(df["decimalLongitude"], df["decimalLatitude"])] , crs = 4326 )
    return gdf

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

    
