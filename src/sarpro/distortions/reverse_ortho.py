# src/sarpro/distortions/reverse_ortho.py

from sarpro.image_processing.grids import generate_pixel_grid
from sarpro.geometry.rpc import inverse_rpc
from sarpro.geometry.transforms import (
    transform_coordinates_between_crs,
    apply_affine_transform
)
from sarpro.utils import get_interpolation_backend
import rasterio.crs

def process_reverse_ortho_tile(
    r_start, r_end, c_start, c_end,
    width, height, z_val, rpc_meta,
    ortho_src, ortho_gpu_data,
    dem_gpu_data=None, dem_inv_transform=None,
    array_module=None
):
    """
    Process reverse orthorectification tile.
    
    This is a WORKFLOW COMPOSITION of:
    - Grid generation
    - Inverse RPC
    - CRS transformation
    - Resampling
    """
    if array_module is None:
        from sarpro.utils import get_array_backend
        array_module, _ = get_array_backend()
    
    cp = array_module
    
    # 1. Generate target grid (ATOMIC)
    grid_r, grid_c = generate_pixel_grid(r_start, r_end, c_start, c_end, array_module=cp)
    
    # 2. Inverse RPC (COMPOSED)
    lon_gpu, lat_gpu = inverse_rpc(
        grid_r, grid_c, z_val, rpc_meta,
        dem_data=dem_gpu_data,
        dem_inv_transform=dem_inv_transform,
        array_module=cp
    )
    
    # 3. Transfer to CPU for CRS transform
    lon_cpu = cp.asnumpy(lon_gpu)
    lat_cpu = cp.asnumpy(lat_gpu)
    
    # 4. Transform CRS (ATOMIC)
    wgs84_crs = rasterio.crs.CRS.from_epsg(4326)
    xs, ys = transform_coordinates_between_crs(
        lon_cpu, lat_cpu, wgs84_crs, ortho_src.crs
    )
    
    # 5. Apply affine transform (ATOMIC)
    u_cpu, v_cpu = apply_affine_transform(xs, ys, ortho_src.transform)
    
    # 6. Sample ortho image (ATOMIC)
    u_gpu = cp.asarray(u_cpu, dtype=cp.float32)
    v_gpu = cp.asarray(v_cpu, dtype=cp.float32)
    sample_coords = cp.stack([v_gpu, u_gpu])
    
    map_coordinates, _ = get_interpolation_backend()
    
    output_bands = []
    for band_idx in range(ortho_gpu_data.shape[0]):
        sampled = map_coordinates(
            ortho_gpu_data[band_idx],
            sample_coords,
            order=1,
            mode='constant',
            cval=0
        )
        output_bands.append(sampled)
    
    return cp.stack(output_bands)
