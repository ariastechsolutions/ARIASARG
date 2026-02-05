from __future__ import annotations

from typing import Optional, Tuple, List, Dict

import rasterio

from sarpro.geometry.transforms import transform_bbox_to_wgs84
from sarpro.io.main import get_tile_info
from sarpro.dem.core import (
    download_clipped_dem_single_tile,
    download_clipped_dem_multi_tile,
    download_full_dem_tiles_no_merge,
    download_full_dem_tiles_merge,
)


def download_copernicus_dem(
    output_path: str,
    *,
    mode: int = 0,
    bbox: Optional[tuple[float, float, float, float]] = None,
    source_crs: str = "EPSG:4326",
    image_path: Optional[str] = None,
    nodata_value: float = -9999.0,
) -> Dict:
    """
    High-level DEM downloader for Copernicus DEM 10 m.

    This function handles:
    - Getting the AOI from either:
        * a raster image (image extent + CRS), OR
        * a user-provided bbox (+ CRS)
    - Transforming AOI to WGS84
    - Getting Copernicus tile list covering the AOI
    - Choosing the correct download strategy based on:
        * mode (0 = clipped, 1 = full tiles)
        * number of tiles (single vs multiple)

    Parameters
    ----------
    output_path : str
        Path where the DEM GeoTIFF will be saved.
    mode : int, optional
        0 = clip to AOI bbox
        1 = full tile(s) without clipping
    bbox : tuple, optional
        (min_x, min_y, max_x, max_y) in `source_crs`.
        If `source_crs == "EPSG:4326"`, this is already WGS84.
        Ignored if `image_path` is provided.
    source_crs : str, optional
        CRS of the provided bbox (e.g., "EPSG:32629").
        Defaults to "EPSG:4326".
    image_path : str, optional
        Path to a raster image. Its extent + CRS are used as AOI.
        If provided, this takes precedence over `bbox`.
    nodata_value : float, optional
        Nodata value for clipped DEMs (mode 0).

    Returns
    -------
    dict
        A dictionary containing:
        - "mode"          : int
        - "bbox_wgs84"    : (min_lon, min_lat, max_lon, max_lat)
        - "tiles"         : list of tile names
        - plus any stats returned by the underlying download functions
    """
    # ------------------------------------------------------------
    # 1. Determine AOI bbox in WGS84
    # ------------------------------------------------------------
    if image_path is not None:
        # Use image extent + CRS
        with rasterio.open(image_path) as src:
            bounds = src.bounds  # left, bottom, right, top
            src_crs_str = src.crs.to_string()

        bbox_src = (bounds.left, bounds.bottom, bounds.right, bounds.top)
        bbox_wgs84 = transform_bbox_to_wgs84(bbox_src, src_crs_str)
        print(f"AOI from image extent in {src_crs_str}: {bbox_src}")
    else:
        if bbox is None:
            raise ValueError("Either 'image_path' or 'bbox' must be provided.")

        # Use provided bbox + CRS
        if source_crs == "EPSG:4326":
            bbox_wgs84 = bbox
        else:
            bbox_wgs84 = transform_bbox_to_wgs84(bbox, source_crs)
        print(f"AOI from bbox in {source_crs}: {bbox}")

    min_lon, min_lat, max_lon, max_lat = bbox_wgs84
    print(
        f"AOI in WGS84: "
        f"[{min_lon:.6f}, {min_lat:.6f}, {max_lon:.6f}, {max_lat:.6f}]"
    )

    # ------------------------------------------------------------
    # 2. Get tile list covering AOI
    # ------------------------------------------------------------
    tiles = get_tile_info(bbox_wgs84)
    if not tiles:
        raise RuntimeError("No Copernicus tiles found for the given AOI.")

    print(f"Found {len(tiles)} tile(s) covering the area:")
    for tile_name, _ in tiles:
        print(f"  - {tile_name}")

    # ------------------------------------------------------------
    # 3. Select strategy based on mode + number of tiles
    # ------------------------------------------------------------
    if mode not in (0, 1):
        raise ValueError("mode must be 0 (clipped) or 1 (full tiles).")

    if mode == 0:
        # --- MODE 0: CLIPPED TO AOI ---
        print("Mode 0 – Downloading clipped DEM...")
        if len(tiles) == 1:
            stats = download_clipped_dem_single_tile(
                tile_urls=tiles,
                bbox=bbox_wgs84,
                output_path=output_path,
                nodata_value=nodata_value,
            )
        else:
            stats = download_clipped_dem_multi_tile(
                tile_urls=tiles,
                bbox=bbox_wgs84,
                output_path=output_path,
                nodata_value=nodata_value,
            )
    else:
        # --- MODE 1: FULL TILE(S) ---
        print("Mode 1 – Downloading full tile(s) (no clipping)...")
        if len(tiles) == 1:
            stats = download_full_dem_tiles_no_merge(
                tile_urls=tiles,
                output_path=output_path,
            )
        else:
            stats = download_full_dem_tiles_merge(
                tile_urls=tiles,
                output_path=output_path,
            )

    # ------------------------------------------------------------
    # 4. Return metadata + stats
    # ------------------------------------------------------------
    result: Dict = {
        "mode": mode,
        "bbox_wgs84": bbox_wgs84,
        "tiles": [t[0] for t in tiles],
    }
    if isinstance(stats, dict):
        result.update(stats)

    print(f"DEM saved to: {output_path}")
    return result
