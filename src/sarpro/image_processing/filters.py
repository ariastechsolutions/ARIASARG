try:
    import cupy as cp
    from cupyx.scipy.ndimage import uniform_filter
    print("Using GPU acceleration with CuPy.")
except ImportError:
    import numpy as cp
    from scipy.ndimage import uniform_filter



def lee_filter_gpu(img, win_size=3, enl=1.0):
    """
    Applies a Lee Speckle Filter to a SAR image using GPU acceleration.

    Args:
        img (numpy.ndarray or cupy.ndarray): Input SAR image (Intensity or Amplitude).
        win_size (int): The size of the window (e.g., 3, 5, 7).
        enl (float): Equivalent Number of Looks. 
                     For a standard intensity image, this is often 1.0.
                     For multi-looked images, this value should be higher.

    Returns:
        numpy.ndarray: The filtered image (on CPU).
    """
    # 1. Ensure input is on the GPU
    # If it's already a cupy array, this is a no-op.
    img_gpu = cp.asarray(img, dtype=cp.float32)

    # 2. Define Noise Variance parameters
    # For intensity images, the coefficient of variation of noise (Cu) 
    # is related to ENL by: Cu^2 = 1 / ENL
    cu_squared = 1.0 / enl

    # 3. Calculate Local Statistics (Mean and Mean of Squares)
    # We use uniform_filter to get the local mean over the window
    local_mean = uniform_filter(img_gpu, size=win_size)
    
    # Calculate local mean of the square of the image (E[x^2])
    local_sqr_mean = uniform_filter(img_gpu**2, size=win_size)

    # 4. Calculate Local Variance
    # Var(x) = E[x^2] - (E[x])^2
    local_var = local_sqr_mean - local_mean**2

    # 5. Calculate the Lee Weighting Factor (K)
    # Formula: K = 1 - (Cu^2 / Ci^2)
    # Where Ci is the image coefficient of variation: Ci = sigma / mean
    # Rearranged for stability: K = 1 - ( (Cu^2 * mean^2) / var )
    
    # We use a small epsilon to avoid division by zero in perfectly flat areas
    epsilon = 1e-10
    
    # Calculate the theoretical variance due to speckle alone
    # speckle_var = (mean^2) * Cu^2
    speckle_var = (local_mean**2) * cu_squared
    
    # Calculate weights
    # We clip the weight to [0, 1] to handle cases where local variance < speckle variance
    # (which happens in homogeneous areas due to estimation error)
    weights = 1.0 - (speckle_var / (local_var + epsilon))
    weights = cp.clip(weights, 0.0, 1.0)

    # 6. Apply the Filter
    # Pixel_out = Mean + K * (Pixel_in - Mean)
    filtered_gpu = local_mean + weights * (img_gpu - local_mean)

    # 7. Return to CPU (Numpy)
    return cp.asnumpy(filtered_gpu)
