import numpy as np
from typing import Union, List
from sarpro.utils import get_array_backend

def precompute_polynomial_powers(
    P: Union[float, np.ndarray],
    L: Union[float, np.ndarray],
    H: Union[float, np.ndarray],
) -> dict:
    """
    Precompute powers of normalized coordinates for RPC polynomial.
    
    Returns
    -------
    dict
        Dictionary with keys: P, L, H, P2, L2, H2, P3, L3, H3
    """
    
    L2 = L * L
    L3 = L2 * L
    P2 = P * P
    P3 = P2 * P
    H2 = H * H
    H3 = H2 * H
    
    return {
        'P': P, 'L': L, 'H': H,
        'P2': P2, 'L2': L2, 'H2': H2,
        'P3': P3, 'L3': L3, 'H3': H3
    }


def evaluate_rpc_polynomial(
    coeffs: List[float],
    powers: dict,
):
    """
    Evaluate a single RPC polynomial using precomputed powers.
    
    Standard 20-term RPC polynomial:
    1, L, P, H, LP, LH, PH, L², P², H², PLH, L³, LP², LH², L²P, P³, PH², L²H, P²H, H³
    
    Parameters
    ----------
    coeffs : list
        20 polynomial coefficients
    powers : dict
        Precomputed powers from precompute_polynomial_powers()
    
    Returns
    -------
    array_like
        Polynomial value
    """
    
    cp, use_gpu = get_array_backend()
    
    c = cp.asarray(coeffs, dtype=cp.float32)
    P, L, H = powers['P'], powers['L'], powers['H']
    P2, L2, H2 = powers['P2'], powers['L2'], powers['H2']
    P3, L3, H3 = powers['P3'], powers['L3'], powers['H3']
    
    val = (c[0]       + c[1]*L     + c[2]*P     + c[3]*H     +
           c[4]*L*P   + c[5]*L*H   + c[6]*P*H   + c[7]*L2    +
           c[8]*P2    + c[9]*H2    + c[10]*P*L*H + c[11]*L3   +
           c[12]*L*P2 + c[13]*L*H2 + c[14]*L2*P + c[15]*P3   +
           c[16]*P*H2 + c[17]*L2*H + c[18]*P2*H + c[19]*H3)
    
    return val


def evaluate_rpc_rational_function(
    num_coeffs: List[float],
    den_coeffs: List[float],
    P,
    L,
    H,
    
):
    """
    Evaluate RPC rational function: numerator / denominator.
    
    Returns
    -------
    array_like
        Normalized coordinate (row_n or col_n)
    """
    powers = precompute_polynomial_powers(P, L, H)
    
    num = evaluate_rpc_polynomial(num_coeffs, powers)
    den = evaluate_rpc_polynomial(den_coeffs, powers)
    
    return num / den