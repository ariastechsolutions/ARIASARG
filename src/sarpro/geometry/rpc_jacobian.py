# src/sarpro/geometry/rpc_jacobian.py

import numpy as np
from typing import Callable, Tuple
from rpc_utils import denormalize_coordinates
from sarpro.utils import get_array_backend

def compute_jacobian_finite_difference(
    P,
    L,
    curr_height,
    rpc_meta: dict,
    forward_func: Callable,
    delta: float = 1e-4,
    array_module=None
) :
    """
    Compute Jacobian matrix using finite differences.
    
    Computes partial derivatives:
        dr/dP, dc/dP (how row/col change with latitude)
        dr/dL, dc/dL (how row/col change with longitude)
    
    Parameters
    ----------
    P, L : array_like
        Current normalized lat/lon estimates
    curr_height : array_like
        Current height values
    rpc_meta : dict
        RPC metadata
    forward_func : callable
        Function signature: forward_func(lon, lat, height, rpc_meta) -> (row, col)
    delta : float
        Perturbation for finite difference
    
    Returns
    -------
    tuple
        (dr_dP, dc_dP, dr_dL, dc_dL) - Jacobian components
    """
    
    
    lat_off = rpc_meta['lat_off']
    lat_scale = rpc_meta['lat_scale']
    lon_off = rpc_meta['long_off']
    lon_scale = rpc_meta['long_scale']
    
    # Current estimates in degrees
    curr_lat = P * lat_scale + lat_off
    curr_lon = L * lon_scale + lon_off
    
    # Forward projection at current point
    r_est, c_est = forward_func(curr_lon, curr_lat, curr_height, rpc_meta)
    
    # Perturb latitude
    lat_p_delta = (P + delta) * lat_scale + lat_off
    r_dp, c_dp = forward_func(curr_lon, lat_p_delta, curr_height, rpc_meta)
    dr_dP = (r_dp - r_est) / delta
    dc_dP = (c_dp - c_est) / delta
    
    # Perturb longitude
    lon_l_delta = (L + delta) * lon_scale + lon_off
    r_dl, c_dl = forward_func(lon_l_delta, curr_lat, curr_height, rpc_meta)
    dr_dL = (r_dl - r_est) / delta
    dc_dL = (c_dl - c_est) / delta
    
    return dr_dP, dc_dP, dr_dL, dc_dL


def solve_newton_step(
    r_err,
    c_err,
    dr_dP,
    dc_dP,
    dr_dL,
    dc_dL,
    array_module=None
):
    """
    Solve Newton-Raphson step using Cramer's rule.
    
    Solves: J * [dP, dL]ᵀ = [r_err, c_err]ᵀ
    
    Returns
    -------
    tuple
        (dP, dL) - Updates to normalized coordinates
    """
    cp , use_gpu = get_array_backend()
    
    # Compute determinant
    det = (dr_dP * dc_dL) - (dr_dL * dc_dP)
    
    # Avoid division by zero
    det = cp.where(cp.abs(det) < 1e-10, 1e-10, det)
    
    # Solve using Cramer's rule
    dP = (dc_dL * r_err - dr_dL * c_err) / det
    dL = (dr_dP * c_err - dc_dP * r_err) / det
    
    return dP, dL
