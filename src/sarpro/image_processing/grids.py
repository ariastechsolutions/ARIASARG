# src/sarpro/image_processing/grids.py

import numpy as np
from typing import Tuple
from sarpro.utils import get_array_backend

def generate_pixel_grid(
    r_start: int,
    r_end: int,
    c_start: int,
    c_end: int,
    pixel_center_offset: float = 0.5,
    array_module=None
) :
    """
    Generate meshgrid of pixel coordinates.
    
    Parameters
    ----------
    r_start, r_end : int
        Row range
    c_start, c_end : int
        Column range
    pixel_center_offset : float
        Offset to sample pixel centers (default: 0.5)
    
    Returns
    -------
    tuple
        (grid_r, grid_c) - 2D arrays of row/col coordinates
    """
    cp, use_gpu = get_array_backend(array_module)
    
    r_coords = cp.arange(r_start, r_end, dtype=cp.float32) + pixel_center_offset
    c_coords = cp.arange(c_start, c_end, dtype=cp.float32) + pixel_center_offset
    
    grid_r, grid_c = cp.meshgrid(r_coords, c_coords, indexing='ij')
    
    return grid_r, grid_c
