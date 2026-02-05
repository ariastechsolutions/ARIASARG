import numpy as np
from typing import Union
from sarpro.utils import (get_array_backend ,
                           get_interpolation_backend)




def normalize_core(
        co_ord_1 : float ,
        co_ord_off : float ,
        co_ord_scale : float
):
    return (co_ord_1 - co_ord_off) / co_ord_scale


def normalize_coordinates(
    lon,
    lat,
    height,
    rpc_meta: dict,
) :
    """
    Normalize geographic coordinates using RPC offsets and scales.
    
    Converts lon/lat/height to normalized coordinates P, L, H:
        L = (lon - lon_off) / lon_scale
        P = (lat - lat_off) / lat_scale
        H = (height - height_off) / height_scale
    
    Returns
    -------
    tuple
        (P, L, H) - Normalized coordinates
    """
    
    lon_off = rpc_meta['long_off']
    lon_scale = rpc_meta['long_scale']
    lat_off = rpc_meta['lat_off']
    lat_scale = rpc_meta['lat_scale']
    h_off = rpc_meta['height_off']
    h_scale = rpc_meta['height_scale']
    
    
    P = normalize_core(lat , lat_off, lat_scale)
    L = normalize_core(lon , lon_off, lon_scale)
    H = normalize_core(height , h_off, h_scale)
    
    return P, L, H






def denormalize_coordinates(
    P,
    L,
    rpc_meta: dict,
    
):
    """
    Denormalize P, L to geographic coordinates.
    
    Returns
    -------
    tuple
        (lon, lat) in degrees
    """
 

    lat_off = rpc_meta['lat_off']
    lat_scale = rpc_meta['lat_scale']
    lon_off = rpc_meta['long_off']
    lon_scale = rpc_meta['long_scale']
    
    lat = P * lat_scale + lat_off
    lon = L * lon_scale + lon_off
    
    return lon, lat


def denormalize_image_coordinates(
    row_n,
    col_n,
    rpc_meta: dict,
):
    """
    Denormalize image coordinates from RPC normalized space to pixels.
    
    Returns
    -------
    tuple
        (row, col) in pixel coordinates
    """
    
    line_off = rpc_meta['line_off']
    line_scale = rpc_meta['line_scale']
    samp_off = rpc_meta['samp_off']
    samp_scale = rpc_meta['samp_scale']
    
    row = row_n * line_scale + line_off
    col = col_n * samp_scale + samp_off
    
    return row, col


