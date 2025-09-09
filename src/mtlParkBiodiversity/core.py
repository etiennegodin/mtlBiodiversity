import geopandas as gpd
from shapely.geometry import Polygon, MultiPolygon


target_crs = 4326

def convert_crs(gdf : gpd.GeoDataFrame, target_crs = 4326, verbose = False):
    """
    Converts gdf to taregt crs 
    """
    if gdf.crs != target_crs:
        try:
            gdf = gdf.to_crs(target_crs)
        except Exception as e:
            print(f"Failed to convert gdf to taregt crs {e}")

    else:
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

    
