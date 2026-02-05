import rasterio
from rasterio.vrt import WarpedVRT
from rasterio.enums import Resampling
from contextlib import contextmanager
from typing import Callable, Generator, Tuple


def get_array_backend():
    """
    Get array computation backend (CuPy for GPU, NumPy for CPU).
    
    Returns
    -------
    module
        CuPy module if available, otherwise NumPy
    bool
        True if using GPU (CuPy), False if using CPU (NumPy)
    
    Example
    -------
    >>> cp, use_gpu = get_array_backend()
    >>> if use_gpu:
    ...     print("Using GPU acceleration")
    >>> data = cp.array([1, 2, 3])
    """
    try:
        import cupy as cp
        print("Using GPU acceleration with CuPy.")
        return cp, True
        
    except ImportError:
        import numpy as cp
        print("Using CPU acceleration with NumPy.")
        return cp, False
    
    


def get_interpolation_backend():
    """Get interpolation backend (cupyx.scipy or scipy)"""
    try:
        from cupyx.scipy.ndimage import map_coordinates
        import cupy as cp
        print("Using GPU acceleration with CuPy.")
        return map_coordinates, cp
    except ImportError:
        from scipy.ndimage import map_coordinates
        import numpy as cp
        print("Using CPU acceleration with NumPy.")
        return map_coordinates, cp


@contextmanager
def align_raster_to_reference(
    source_path: str,
    reference_dataset: rasterio.DatasetReader,
    resampling: Resampling = Resampling.bilinear,
    nodata: float = 0
):
    """
    Align a raster to match another raster's geometry using WarpedVRT.
    
    This creates a virtual warped raster that matches the reference raster's
    CRS, transform, and dimensions. Useful for ensuring pixel-perfect alignment
    between multiple rasters.
    
    Parameters
    ----------
    source_path : str
        Path to raster to be aligned
    reference_dataset : rasterio.DatasetReader
        Opened reference raster (defines target geometry)
    resampling : Resampling, optional
        Resampling method (default: bilinear)
    nodata : float, optional
        NoData value (default: 0)
    
    Yields
    ------
    rasterio.DatasetReader
        Warped raster matching reference geometry
    
    Example
    -------
    >>> with rasterio.open('reference.tif') as ref:
    ...     with align_raster_to_reference('source.tif', ref) as aligned:
    ...         data = aligned.read(1)
    """
    with rasterio.open(source_path) as src:
        vrt_options = {
            'resampling': resampling,
            'crs': reference_dataset.crs,
            'transform': reference_dataset.transform,
            'height': reference_dataset.height,
            'width': reference_dataset.width,
            'nodata': nodata
        }
        
        with WarpedVRT(src, **vrt_options) as vrt:
            yield vrt



def generate_blocks(
    width: int,
    height: int,
    block_size: int = 4096
) -> Generator[Tuple[rasterio.windows.Window, int, int], None, None]:
    """
    Generate windows for block-wise raster processing.
    
    Parameters
    ----------
    width : int
        Raster width in pixels
    height : int
        Raster height in pixels
    block_size : int, optional
        Block size in pixels (default: 4096)
    
    Yields
    ------
    tuple
        (window, block_number, total_blocks)
    
    Example
    -------
    >>> for window, block_num, total in generate_blocks(10000, 8000, 2048):
    ...     print(f"Processing block {block_num}/{total}")
    """
    total_blocks = ((height + block_size - 1) // block_size) * \
                   ((width + block_size - 1) // block_size)
    
    block_num = 0
    for r in range(0, height, block_size):
        for c in range(0, width, block_size):
            window = rasterio.windows.Window(
                c, r,
                min(block_size, width - c),
                min(block_size, height - r)
            )
            block_num += 1
            yield window, block_num, total_blocks