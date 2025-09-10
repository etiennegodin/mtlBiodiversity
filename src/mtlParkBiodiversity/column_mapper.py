


def unify_columns(gdf, mapping):
    """
    Unify column names to a standard set.
    """
  
    gdf = gdf.rename(columns=column_mapping)
    
    # Ensure all required columns are present
    required_columns = ['name', 'area_m2', 'type', 'geometry']
    for col in required_columns:
        if col not in gdf.columns:
            gdf[col] = None  # or some default value
    
    return gdf[required_columns]