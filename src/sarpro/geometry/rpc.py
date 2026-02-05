# src/sarpro/geometry/rpc.py

from sarpro.geometry.rpc_utils import normalize_coordinates, denormalize_image_coordinates , denormalize_coordinates
from sarpro.geometry.rpc_poynomial import evaluate_rpc_rational_function
from sarpro.geometry.rpc_jacobian import compute_jacobian_finite_difference, solve_newton_step
from sarpro.dem.sampling import sample_dem_at_coordinates

def apply_rpc(lon, lat, height, rpc_meta, array_module=None):
    """
    Apply RPC: ground → image coordinates.
    
    """
    # Step 1: Normalize inputs
    P, L, H = normalize_coordinates(lon, lat, height, rpc_meta)
    
    # Step 2: Evaluate rational functions
    row_n = evaluate_rpc_rational_function(
        rpc_meta['line_num_coeff'],
        rpc_meta['line_den_coeff'],
        P, L, H
    )
    
    col_n = evaluate_rpc_rational_function(
        rpc_meta['samp_num_coeff'],
        rpc_meta['samp_den_coeff'],
        P, L, H,
    )
    
    # Step 3: Denormalize to pixel coordinates
    row, col = denormalize_image_coordinates(row_n, col_n, rpc_meta)
    
    return row, col



def inverse_rpc(
    target_rows, target_cols, fixed_height, rpc_meta,
    dem_data=None, dem_inv_transform=None, iterations=5, array_module=None
):
    """
    Inverse RPC: image → ground coordinates.
    
    This is a COMPOSITION of:
    - Iterative Newton-Raphson solver
    - Jacobian computation
    - Optional DEM sampling
    """
    if array_module is None:
        from sarpro.utils import get_array_backend
        array_module, _ = get_array_backend()
    
    cp = array_module
    
    # Initialize normalized coordinates
    P = cp.zeros_like(target_rows, dtype=cp.float32)
    L = cp.zeros_like(target_cols, dtype=cp.float32)
    
    # Initialize height
    if dem_data is not None:
        curr_height = cp.full_like(target_rows, rpc_meta['height_off'], dtype=cp.float32)
    else:
        curr_height = fixed_height
    
    # Newton-Raphson iterations
    for i in range(iterations):
        # Denormalize current estimate
        curr_lon, curr_lat = denormalize_coordinates(P, L, rpc_meta)
        
        # Update height from DEM if available
        if dem_data is not None and dem_inv_transform is not None:
            curr_height = sample_dem_at_coordinates(
                curr_lon, curr_lat, dem_data, dem_inv_transform, cp
            )
        
        # Forward projection
        r_est, c_est = apply_rpc(curr_lon, curr_lat, curr_height, rpc_meta)
        
        # Residuals
        r_err = target_rows - r_est
        c_err = target_cols - c_est
        
        # Compute Jacobian
        dr_dP, dc_dP, dr_dL, dc_dL = compute_jacobian_finite_difference(
            P, L, curr_height, rpc_meta, apply_rpc, array_module=cp
        )
        
        # Solve Newton step
        dP, dL = solve_newton_step(r_err, c_err, dr_dP, dc_dP, dr_dL, dc_dL, cp)
        
        # Update
        P += dP
        L += dL
    
    # Final denormalization
    final_lon, final_lat = denormalize_coordinates(P, L, rpc_meta, cp)
    
    return final_lon, final_lat