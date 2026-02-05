import rasterio
import numpy as np
from rasterio.merge import merge
from typing import List, Tuple, Optional, Dict


TileList = List[Tuple[str, str]]  # (tile_name, s3_url)


def open_dem_sources(tile_urls: TileList) -> List[rasterio.io.DatasetReader]:
    """Open all DEM tiles from (tile_name, s3_url) list."""
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

    return src_files


def merge_dem_sources(
    src_files: List[rasterio.io.DatasetReader],
    bounds: Optional[tuple[float, float, float, float]] = None,
):
    """
    Merge DEM sources. If bounds is provided, clip to bounds.
    Returns (mosaic, transform).
    """
    if bounds is not None:
        mosaic, transform = merge(src_files, bounds=bounds)
    else:
        mosaic, transform = merge(src_files)

    return mosaic, transform


def compute_dem_stats(data: np.ndarray, valid_mask: np.ndarray) -> Dict[str, float]:
    """
    Compute simple statistics for DEM: valid pixels, total pixels, min/max elevation.
    """
    valid_data = data[valid_mask]
    total_pixels = int(data.size)
    valid_pixels = int(valid_data.size)

    stats: Dict[str, float] = {
        "total_pixels": total_pixels,
        "valid_pixels": valid_pixels,
    }

    if valid_pixels > 0:
        stats["min_elevation"] = float(valid_data.min())
        stats["max_elevation"] = float(valid_data.max())

    return stats


def print_dem_stats(stats: Dict[str, float]) -> None:
    """Pretty-print DEM statistics."""
    valid_pixels = stats["valid_pixels"]
    total_pixels = stats["total_pixels"]

    if valid_pixels > 0:
        print(f"✓ Valid pixels: {valid_pixels:,} / {total_pixels:,}")
        print(
            f"✓ Elevation range: {stats['min_elevation']:.1f} "
            f"to {stats['max_elevation']:.1f} meters"
        )
    else:
        print("⚠ Warning: No valid elevation data in this area")


def write_single_band_dem(
    output_path: str,
    template_profile: dict,
    data: np.ndarray,
    transform,
    nodata_value: Optional[float] = None,
) -> None:
    """Write a single-band DEM with LZW compression."""
    profile = template_profile.copy()
    profile.update(
        {
            "height": data.shape[0],
            "width": data.shape[1],
            "transform": transform,
            "compress": "lzw",
            "dtype": rasterio.float32,
        }
    )

    if nodata_value is not None:
        profile["nodata"] = nodata_value

    with rasterio.open(output_path, "w", **profile) as dst:
        dst.write(data, 1)


def write_multiband_dem(
    output_path: str,
    template_profile: dict,
    data: np.ndarray,
    transform,
) -> None:
    """Write multi-band DEM (e.g. merged full tiles)."""
    profile = template_profile.copy()
    profile.update(
        {
            "height": data.shape[1],
            "width": data.shape[2],
            "transform": transform,
            "compress": "lzw",
        }
    )

    with rasterio.open(output_path, "w", **profile) as dst:
        dst.write(data)
