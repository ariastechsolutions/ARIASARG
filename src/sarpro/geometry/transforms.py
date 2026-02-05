from pyproj import CRS, Transformer
import numpy as np
from pyproj import CRS
from rasterio.warp import transform as rasterio_transform
from typing import Tuple
import rasterio.crs


def transform_bbox(
    bbox: tuple[float, float, float, float],
    source_crs: str ,
    target_crs: str = 'EPSG:4326'
) -> tuple[float, float, float, float]:
    min_x, min_y, max_x, max_y = bbox

    src_crs = CRS.from_user_input(source_crs)
    dst_crs = CRS.from_user_input(target_crs)

    if src_crs == dst_crs:
        return min_x, min_y, max_x, max_y

    transformer = Transformer.from_crs(src_crs, dst_crs, always_xy=True)

    # Transform all four corners
    pts_src = [
        (min_x, min_y),
        (min_x, max_y),
        (max_x, min_y),
        (max_x, max_y),
    ]

    lons, lats = zip(*[transformer.transform(x, y) for x, y in pts_src])

    min_lon = min(lons)
    max_lon = max(lons)
    min_lat = min(lats)
    max_lat = max(lats)

    return (min_lon, min_lat, max_lon, max_lat)



def calculate_geometric_factor(eia_deg, lia_rad, array_module=None):
    """
    Compute RTC geometric normalization factor.
    
    The geometric factor is the ratio of local to ellipsoid incidence angles:
        gamma = sin(theta_local) / sin(theta_ellipsoid)
    
    This factor corrects for terrain-induced radiometric distortions in SAR imagery.
    
    Parameters
    ----------
    eia_deg : array_like
        Ellipsoid Incidence Angle in DEGREES
    lia_rad : array_like
        Local Incidence Angle in RADIANS
    array_module : module, optional
        Computation backend (cupy or numpy). If None, auto-detected.
    
    Returns
    -------
    array_like
        Geometric factor (same shape as inputs), values >= 0
    
    Notes
    -----
    - Division by zero is avoided by clamping sin(theta_ell) to minimum 1e-6
    - Negative values are clamped to 0.0
    - Works with both NumPy (CPU) and CuPy (GPU) arrays
    
    References
    ----------
    Small, D. (2011). Flattening Gamma: Radiometric Terrain Correction for SAR Imagery.
    IEEE Transactions on Geoscience and Remote Sensing, 49(8), 3081-3093.
    """
    if array_module is None:
        from sarpro.utils.general import get_array_backend
        array_module, _ = get_array_backend()
    
    cp = array_module
    
    # Convert EIA from degrees to radians
    eia_rad = cp.radians(eia_deg)
    
    # Calculate trigonometric terms
    sin_theta_loc = cp.sin(lia_rad)
    sin_theta_ell = cp.sin(eia_rad)
    
    # Avoid division by zero
    sin_theta_ell = cp.where(sin_theta_ell < 1e-6, 1e-6, sin_theta_ell)
    
    # Calculate geometric factor
    ratio = sin_theta_loc / sin_theta_ell
    
    # Clamp negative values to zero
    ratio = cp.maximum(ratio, 0.0)
    
    return ratio


def apply_rtc_correction(
    dn, 
    geometric_factor, 
    k_factor, 
    output_db=False,
    array_module=None
):
    """
    Apply Radiometric Terrain Correction to SAR Digital Numbers.
    
    Converts raw SAR DN values to radiometrically terrain-corrected
    backscatter (sigma0).
    
    Formula:
        sigma_raw = DN^2 / K
        sigma_rtc = sigma_raw * geometric_factor
    
    Parameters
    ----------
    dn : array_like
        SAR Digital Number values
    geometric_factor : array_like
        Pre-computed geometric factor (from calculate_geometric_factor)
    k_factor : float
        Calibration constant (sensor-specific)
    output_db : bool, optional
        If True, output in dB scale (10*log10). Default: False (linear)
    array_module : module, optional
        Computation backend (cupy or numpy)
    
    Returns
    -------
    array_like
        RTC-corrected backscatter values
    
    Notes
    -----
    - Negative values are clamped to 0.0
    - For dB output, zero/negative values are set to -99.0 dB
    """
    if array_module is None:
        from utils.general import get_array_backend
        array_module, _ = get_array_backend()
    
    cp = array_module
    
    dn = dn.astype(cp.float32)
    
    # Convert DN to raw sigma0
    sigma_raw = (dn ** 2) / k_factor
    
    # Apply geometric correction
    sigma_rtc = sigma_raw * geometric_factor
    
    # Clamp negative values
    sigma_rtc = cp.maximum(sigma_rtc, 0.0)
    
    if output_db:
        # Convert to dB, avoiding log(0)
        sigma_rtc = cp.where(sigma_rtc <= 0, -99.0, 10 * cp.log10(sigma_rtc))
    
    return sigma_rtc


def ecef_to_enu_rotation_matrix(lat_deg, lon_deg , array_module=None):
    """
    Constructs the rotation matrix to project ECEF vectors into 
    the Local Tangent Plane (East-North-Up) centered at a specific lat/lon.
    """
    if array_module is None:
        from utils.general import get_array_backend
        array_module, _ = get_array_backend()

    cp = array_module
    rad_lat = cp.radians(lat_deg)
    rad_lon = cp.radians(lon_deg)
    
    sin_lat = cp.sin(rad_lat)
    cos_lat = cp.cos(rad_lat)
    sin_lon = cp.sin(rad_lon)
    cos_lon = cp.cos(rad_lon)

    # Rotation matrix R (ECEF -> ENU)
    # Row 1: East unit vector
    # Row 2: North unit vector
    # Row 3: Up unit vector
    R = cp.array([
        [-sin_lon,             cos_lon,              0],
        [-sin_lat * cos_lon,  -sin_lat * sin_lon,    cos_lat],
        [cos_lat * cos_lon,    cos_lat * sin_lon,    sin_lat]
    ])
    
    return R


def convert_ecef_to_enu(vector_ecef, 
                        scene_center_ecef, 
                        scene_geo,
                        is_position=True,
                        array_module=None):
    """
    Converts a vector from ECEF to ENU coordinates centered at scene_center_ecef.
    
    Args:
        vector_ecef: The vector to convert [X, Y, Z].
        scene_center_ecef: The origin of the ENU frame [X, Y, Z].
        scene_geo: [Lat, Lon] of the scene center (degrees).
        is_position: True if vector_ecef is a position (subtracts center), 
                     False if it's a direction/velocity vector.
        
    Returns:
        np.array: Vector in ENU coordinates [East, North, Up].
    """
    if array_module is None:
        from utils.general import get_array_backend
        array_module, _ = get_array_backend()
    cp = array_module
    # 1. Calculate the Vector to project
    if is_position:
        # For position: vector points FROM the Center TO the point
        vec_to_project = cp.array(vector_ecef) - cp.array(scene_center_ecef)
    else:
        # For velocity/direction: just use the vector as is
        vec_to_project = cp.array(vector_ecef)
    
    # 2. Get Rotation Matrix based on Scene Center Lat/Lon
    # Fix: Ensure we pass lat (index 0) and lon (index 1) explicitly
    R = ecef_to_enu_rotation_matrix(scene_geo[0], scene_geo[1], array_module=cp)
    
    # 3. Project Vector
    vec_enu = R.dot(vec_to_project)
    
    return vec_enu



def enu_to_elevation_angle(east : float,
                            north : float,
                            up : float,
                            array_module=None):
    """
    Converts ENU coordinates (meters) to Elevation Angle (degrees).
    
    Args:
        east (float): East component (meters)
        north (float): North component (meters)
        up (float): Up/Altitude component (meters)
        
    Returns:
        float: Angle of Altitude (Elevation) in degrees.
               90 = Zenith, 0 = Horizon.
    """
    if array_module is None:
        from utils.general import get_array_backend
        array_module, _ = get_array_backend()
    cp = array_module
    # 1. Calculate Horizontal Ground Distance
    ground_dist = cp.hypot(east, north)
    
    # 2. Calculate Angle (Result is in Radians)
    # arctan2(y, x) -> arctan2(vertical, horizontal)
    elevation_rad = cp.arctan2(up, ground_dist)
    
    # 3. Convert to Degrees
    elevation_deg = cp.degrees(elevation_rad)
    
    return elevation_deg


def transform_coordinates_between_crs(
    x: np.ndarray,
    y: np.ndarray,
    source_crs: rasterio.crs.CRS,
    target_crs: rasterio.crs.CRS
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Transform arrays of coordinates between CRS systems.
    
    Low-level function for bulk coordinate transformation.
    """
    if source_crs == target_crs:
        return x, y
    
    shape = x.shape
    x_flat, y_flat = rasterio_transform(
        source_crs, target_crs,
        x.ravel(), y.ravel()
    )
    
    x_out = np.array(x_flat).reshape(shape)
    y_out = np.array(y_flat).reshape(shape)
    
    return x_out, y_out

def apply_affine_transform(
    x,
    y,
    affine_transform
) :
    """
    Apply affine transformation to map coordinates → pixel coordinates.
    
    Parameters
    ----------
    x, y : array_like
        Map coordinates (meters or degrees)
    affine_transform : rasterio.Affine
        Affine transformation matrix
    
    Returns
    -------
    tuple
        (col, row) in pixel coordinates
    """
    inv = ~affine_transform
    
    col = x * inv.a + y * inv.b + inv.c
    row = x * inv.d + y * inv.e + inv.f
    
    return col, row