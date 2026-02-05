# src/sarpro/dem/sampling.py

import numpy as np
from typing import Tuple, Optional

def sample_dem_at_coordinates(
    lon,
    lat,
    dem_data,
    dem_inv_transform: Tuple[float, float, float, float, float, float],
    array_module=None
) :
    """
    Sample DEM elevation at given geographic coordinates.
    
    Parameters
    ----------
    lon, lat : array_like
        Geographic coordinates in degrees
    dem_data : array_like
        DEM elevation array (height, width)
    dem_inv_transform : tuple
        Inverse affine transform (a, b, c, d, e, f) for lon/lat → pixel
    
    Returns
    -------
    array_like
        Sampled elevation values
    """
    
    
    from sarpro.utils import get_interpolation_backend , get_array_backend
    cp, use_gpu = get_array_backend(array_module)
    map_coordinates, _ = get_interpolation_backend()
    
    a, b, c, d, e, f = dem_inv_transform
    
    # Map lon/lat to DEM pixel coordinates
    dem_col = a * lon + b * lat + c
    dem_row = d * lon + e * lat + f
    
    # Stack for map_coordinates (expects [row, col])
    sample_coords = cp.stack([dem_row, dem_col])
    
    # Sample with bilinear interpolation
    elevation = map_coordinates(dem_data, sample_coords, order=1, mode='nearest')
    
    return elevation
