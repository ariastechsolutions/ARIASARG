import numpy as np
from typing import Tuple



def mask_invalid_dem_values(
    data: np.ndarray,
    nodata_value: float = -9999,
    valid_threshold = 0
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Mask invalid DEM values and replace with nodata.
    
    Parameters:
    -----------
    data : np.ndarray
        DEM elevation data
    nodata_value : float
        Value to assign to invalid pixels
    valid_threshold : float
        Pixels > threshold are considered valid
    
    Returns:
    --------
    Tuple[np.ndarray, np.ndarray]
        (masked_data, valid_mask)
    """
    # Create mask: valid pixels are > valid_threshold
    valid_mask = data > valid_threshold

    # Convert to float32 and apply nodata to invalid pixels
    masked_data = data.astype(np.float32)
    masked_data[~valid_mask] = nodata_value

    return masked_data, valid_mask
