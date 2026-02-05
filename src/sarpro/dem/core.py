import rasterio
from rasterio.windows import from_bounds
from sarpro.image_processing.main import mask_invalid_dem_values
import numpy as np
from rasterio.merge import merge
from sarpro.utils.dem_utils import (
    compute_dem_stats,
    print_dem_stats,
    write_single_band_dem,
)


def download_clipped_dem_single_tile(
    tile_urls: list[tuple[str, str]],
    bbox: tuple[float, float, float, float],
    output_path: str,
    nodata_value: float = -9999,
) -> dict:
    min_lon, min_lat, max_lon, max_lat = bbox
    tile_name, s3_url = tile_urls[0]
    print(f"Reading from single tile: {tile_name}")

    try:
        with rasterio.open(s3_url) as src:
            window = from_bounds(min_lon, min_lat, max_lon, max_lat, src.transform)
            data = src.read(1, window=window)
            transform = src.window_transform(window)

            masked_data, valid_mask = mask_invalid_dem_values(
                data, nodata_value, valid_threshold=0.0
            )

            stats = compute_dem_stats(data, valid_mask)
            print_dem_stats(stats)

            write_single_band_dem(
                output_path=output_path,
                template_profile=src.profile,
                data=masked_data,
                transform=transform,
                nodata_value=nodata_value,
            )

            return stats

    except Exception as e:
        raise Exception(f"Failed to download DEM: {str(e)}")


def download_clipped_dem_multi_tile(
    tile_urls: list[tuple[str, str]],
    bbox: tuple[float, float, float, float],
    output_path: str,
    nodata_value: float = -9999
) -> dict:
    """
    Download, merge, and clip multiple Copernicus DEM tiles.
    
    Parameters:
    -----------
    tile_urls : list[tuple[str, str]]
        List of (tile_name, s3_url) tuples
    bbox : tuple
        (min_lon, min_lat, max_lon, max_lat) in WGS84
    output_path : str
        Path to save merged and clipped DEM
    nodata_value : float
        Value for invalid pixels
    
    Returns:
    --------
    dict
        Statistics: valid_pixels, total_pixels, min_elevation, max_elevation
    """
        # Multiple tiles: merge and crop
    print(f"Merging {len(tile_urls)} tiles (AOI spans multiple 1° tiles)...")
        
    min_lon, min_lat, max_lon, max_lat = bbox



    src_files = []
    for tile_name, s3_url in tile_urls:
        try:
            src = rasterio.open(s3_url)
            src_files.append(src)
            print(f"  Opened: {tile_name}")
        except Exception as e:
            print(f"  Warning: Could not open {tile_name}: {e}")
        
    if not src_files:
        raise Exception("No tiles could be opened!")
        
        # Merge with bounds
    mosaic, transform = merge(
        src_files, 
        bounds=(min_lon, min_lat, max_lon, max_lat)
    )
        
        # Get data (first band)
    data = mosaic[0]
    
        # ⭐ MASK OUT INVALID VALUES
    print("Masking invalid pixels (≤0)...")
        
    masked_data , valid_mask = mask_invalid_dem_values(
            data,
            nodata_value , 
            valid_threshold=0.0
    )
    # Calculate statistics
    valid_data = data[valid_mask]
    
    if len(valid_data) > 0:
        print(f"✓ Valid pixels: {len(valid_data):,} / {data.size:,}")
        print(f"✓ Elevation range: {valid_data.min():.1f} to {valid_data.max():.1f} meters")
    else:
        print("⚠ Warning: No valid elevation data in this area")
    
    # Write output
    profile = src_files[0].profile.copy()
    profile.update({
        'height': masked_data.shape[0],
        'width': masked_data.shape[1],
        'transform': transform,
        'compress': 'lzw',
        'dtype': rasterio.float32,
        'nodata': nodata_value
    })
    
    with rasterio.open(output_path, 'w', **profile) as dst:
        dst.write(masked_data[np.newaxis, :, :])
    
    # Close sources
    for src in src_files:
        src.close()
    
    

def download_full_dem_tiles_no_merge(
    tile_urls: list[tuple[str, str]],
    output_path: str
) -> dict:
    """
    Download complete Copernicus DEM tile(s) without clipping.
    
    Parameters:
    -----------
    tile_urls : list[tuple[str, str]]
        List of (tile_name, s3_url) tuples
    output_path : str
        Path to save full tile(s)
    
    Returns:
    --------
    dict
        Dimensions: width, height
    """

    tile_name, s3_url = tile_urls[0]
    print(f"Downloading: {tile_name}")
            
    try:
        with rasterio.open(s3_url) as src:
            data = src.read()
            
            profile = src.profile.copy()
            profile.update({
                'compress': 'lzw'
            })
            
            with rasterio.open(output_path, 'w', **profile) as dst:
                dst.write(data)
            
            print(f"✓ Downloaded full tile: {src.width} × {src.height} pixels")
            
    except Exception as e:
        raise Exception(f"Failed to download tile: {str(e)}")



def download_full_dem_tiles_merge(
    tile_urls: list[tuple[str, str]],
    output_path: str
) -> dict:
    """
    Download complete Copernicus DEM tile(s) and merging without clipping .
    
    Parameters:
    -----------
    tile_urls : list[tuple[str, str]]
        List of (tile_name, s3_url) tuples
    output_path : str
        Path to save full tile(s)
    
    Returns:
    --------
    dict
        Dimensions: width, height
    """


    print(f"Downloading and merging {len(tile_urls)} full tiles...")
    
    src_files = []
    for tile_name, s3_url in tile_urls:
        try:
            src_files.append(rasterio.open(s3_url))
        except Exception as e:
            print(f"  Warning: Could not open {tile_name}: {e}")
    
    if not src_files:
        raise Exception("No tiles could be opened!")
    
    # Merge without bounds restriction
    mosaic, transform = merge(src_files)
    
    profile = src_files[0].profile.copy()
    profile.update({
        'height': mosaic.shape[1],
        'width': mosaic.shape[2],
        'transform': transform,
        'compress': 'lzw'
    })
    
    with rasterio.open(output_path, 'w', **profile) as dst:
        dst.write(mosaic)
    
    for src in src_files:
        src.close()
    
    print(f"✓ Merged full tiles: {mosaic.shape[2]} × {mosaic.shape[1]} pixels")


