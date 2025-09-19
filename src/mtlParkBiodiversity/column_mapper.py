import geopandas as gpd
import json
from pathlib import Path    


def load_column_mapping(mapping_file : Path = None):
    """
    Load column mapping from a JSON file.
    """
    if mapping_file is not None and Path(mapping_file).exists():
        with open(str(mapping_file), 'r') as f:
            mapping = json.load(f)
        return mapping
    return {}

def save_mapping(mapping : json, mapping_file : Path = None):
    """Save mapping to config file."""
    with open(str(mapping_file), "w") as f:
        json.dump(mapping, f, indent=2)

def get_user_mapping(missing_col, available_cols):
    """Ask user which column to map to the missing one."""
    print(f"\nMissing expected column: {missing_col}")
    print("Available columns:", available_cols)
    choice = input(f"Select a column to map to '{missing_col}' (or press Enter to skip): ")
    return choice if choice in available_cols else None

def get_user_expected_col():
    print("Expected columns not provided. Please type wanted columns as list:")
    expected_input = input()
    expected_columns = list(expected_input.split(sep = ","))
    if not isinstance(expected_columns, list):
        raise UserWarning('Could not ')
        get_user_expected_col()
    else:
        return expected_columns


def unify_columns(file : Path, expected_columns : list = None, force = False):
    """
    Unify column names to a standard set.
    """

    print(f'Mapping {file} columns ')
    config_path = file.parent.parent.parent / '.config'

    #Create config directory if not existing
    Path.mkdir(config_path, exist_ok= True)

    # Define mapping file path
    mapping_file = config_path / f"{file.stem}_column_mapping.json"

    # If force, delete mapping to start over
    if force and mapping_file.exists():
        mapping_file.unlink()

    # Load existing mapping or create a new one
    mapping = load_column_mapping(mapping_file = mapping_file)

    # Read the file to get available columns
    gdf = gpd.read_file(file)    
    available_cols = gdf.columns.tolist()

    if expected_columns is None:

            
    # Loop over expected cols
    for col in expected_columns:
        # Check if in gfd
        if col not in available_cols:
            # Check if in mapping first
            if col in mapping:
                #Add to mapping (to keep it while it saves again) 
                mapped_col = mapping[col]

                #Check if previous mapping is still in available columns
                if mapped_col in available_cols:
                    gdf.rename(columns={mapped_col: col}, inplace=True)

                else:
                    print(f"Mapped column '{mapped_col}' for '{col}' not found in available columns.")
            else:
                # Col not in mapping, ask user for input 
                mapped_col = get_user_mapping(col, available_cols)
                if mapped_col:
                    # Remap if provided
                    gdf.rename(columns={mapped_col: col}, inplace=True)
                    mapping[col] = mapped_col
                else:
                    # Not provided by user, skip
                    print(f"No mapping provided for '{col}'. It will be skipped.")


    save_mapping(mapping, mapping_file = mapping_file)

    return gdf