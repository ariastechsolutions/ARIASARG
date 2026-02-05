import math
import os
import ast
from dotenv import dotenv_values


def get_tile_info(bbox: tuple[float, float, float, float]) -> list[tuple[str, str]]:
    """
    Generate Copernicus DEM tile names and S3 URLs for a bounding box.
    
    Parameters:
    -----------
    bbox : tuple[float, float, float, float]
        Bounding box as (min_lon, min_lat, max_lon, max_lat) in WGS84
    
    Returns:
    --------
    list[tuple[str, str]]
        List of (tile_name, s3_url) tuples
    """
    min_lon, min_lat, max_lon, max_lat = bbox
    
    tiles = []
    for lat in range(math.floor(min_lat), math.ceil(max_lat)):
        for lon in range(math.floor(min_lon), math.ceil(max_lon)):
            lat_str = f"N{abs(lat):02d}_00" if lat >= 0 else f"S{abs(lat):02d}_00"
            lon_str = f"E{abs(lon):03d}_00" if lon >= 0 else f"W{abs(lon):03d}_00"
            
            tile_name = f"Copernicus_DSM_COG_10_{lat_str}_{lon_str}_DEM"
            s3_url = f"s3://copernicus-dem-30m/{tile_name}/{tile_name}.tif"
            
            tiles.append((tile_name, s3_url))
    
    return tiles





def load_sar_env(env_path='.env'):
    """
    Loads variables from a .env file into a dictionary using python-dotenv.
    Parses vector strings (e.g., "[1.0 2.0 3.0]") into lists if possible.
    """
    if not os.path.exists(env_path):
        print(f"Warning: {env_path} not found.")
        return {}

    env_vars = dotenv_values(env_path)
    config = {}

    for key, value in env_vars.items():
        if value is None:
            config[key] = None
            continue

        # Try to parse numeric lists (vectors) which look like "[x y z]"
        if value.startswith('[') and value.endswith(']'):
            content = value[1:-1].strip()
            # Handle space-separated items inside brackets
            if ' ' in content and ',' not in content:
                try:
                    config[key] = [float(p) for p in content.split()]
                    continue
                except ValueError:
                    pass
            # Fallback to ast.literal_eval for standard list formats
            try:
                config[key] = ast.literal_eval(value)
                continue
            except (ValueError, SyntaxError):
                pass

        # Try to parse simple floats/ints
        try:
            config[key] = float(value) if '.' in value else int(value)
        except ValueError:
            config[key] = value

    return config

