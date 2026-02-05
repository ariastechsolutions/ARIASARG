# SAR Processing Library (sarpro)

## Introduction

This library provides tools and building blocks for SAR (Synthetic Aperture Radar) processing techniques, including geometry, image processing, and distortion handling. The project is modular so you can implement and document each component independently.


## Table of Contents
- Modules
- Installation
- Quick Start
- Usage Examples
- Contributing
- Roadmap
- License

## Modules

Below are the main modules in the library. 

The Main Entry Point for each module is the Main.py file which is located at the core of each module folder. 

The `main.py` file serves as the primary entry point for each module. It orchestrates the functions and classes implemented across the other files within the module directory.

Accordingly, the documentation focuses on how to use `main.py` for each module and highlights the core functionalities exposed by that module.

### DEM Module
The **DEM Module** handles all tasks related to Digital Elevation Models (DEM), from data acquisition to preprocessing and preparation for integration with SAR products.

#### Current Functionality
- Downloading **Copernicus DEM 30 m** data using:
  - A user-provided bounding box, or
  - A SAR product image, provided that the image contains valid metadata and bounding box information
- Automatic downloading of all DEM tiles required to cover the target area
- Automatic merging of multiple DEM tiles into a single raster
- Full control over:
  - Clipping behavior
  - Tile merging
  - Output format of the final DEM product

#### Core Function

The main functionality is provided through the `download_copernicus_dem` function:

```python
download_copernicus_dem(
    output_path,
    bbox,
    image_path,
    source_crs="EPSG:4326",
    mode
)
```
  #### Parameters

- **output_path**  
  Path where the final DEM raster will be saved.

- **bbox**  
  Bounding box defining the area of interest.  
  If provided, it will be used directly to determine the DEM coverage.

- **image_path**  
  Path to a SAR product image.  
  If provided, the bounding box will be extracted from the image metadata.

- **source_crs** *(default: `EPSG:4326`)*  
  Coordinate Reference System (CRS) to which the DEM will be converted.

- **mode**  
  Defines the DEM download behavior:
  - `"clip"`: Download and clip DEM data to the bounding box
  - `"full"`: Download full DEM tiles covering the bounding box without clipping

Either **bbox** or **image_path** must be provided.  
If multiple DEM tiles are required, merging is performed automatically.
### Geometry
The **Geometry Module** handles geometric transformations, angle calculations, and radiometric terrain corrections (RTC) for SAR data. It provides essential tools for converting raw SAR measurements into radiometrically corrected backscatter values by accounting for terrain-induced distortions.
- **Coordinate System Transformations**: Transform bounding boxes between different coordinate reference systems (CRS) to WGS84 (EPSG:4326)

- **Geometric Factor Computation**: Calculate RTC geometric normalization factors from Ellipsoid Incidence Angle (EIA) and Local Incidence Angle (LIA) data

- **Radiometric Terrain Correction (RTC)**: Apply terrain correction to raw SAR Digital Numbers (DN) to produce radiometrically corrected backscatter values

#### Core Functions:

1 - This function is responsible for computing the RTC geometric factor raster from EIA and LIA rasters. The scientific basis for this computation is detailed in the Knowlegde layer of the Library. 

```python
compute_geometric_factor_raster(
    eia_path: str,
    lia_path: str,
    output_path: str,
    block_size: int ,
    verbose: bool 
)
```
#### Parameters

- **eia_path**  
  Path to the Ellipsoid Incidence Angle (EIA) raster.

- **lia_path**  
  Path to the Local Incidence Angle (LIA) raster.

- **output_path**  
  Path where the final geometric factor raster will be saved.

- **block_size** *(default: `4096`)*  
  Size of processing blocks in pixels of raster for memory management.

- **verbose** *(default: `True`)*  
  If `True`, print verbose output during processing.

 2- This function is responsible for applying the RTC correction to a raw SAR raster using a pre-computed geometric factor raster. The scientific basis for this computation is detailed in the Knowlegde layer of the Library. 

 ```python
apply_rtc_to_raster(
    sar_path: str,
    geometric_factor_path: str,
    k_factor: float,
    output_path: str,
    output_db: bool,
    block_size: int ,
    verbose: bool
)
 ```
#### Parameters

- **sar_path**  
  Path to raw SAR raster (DN or amplitude values)

- **geometric_factor_path**  
  Path to geometric factor raster (can be computed from compute_geometric_factor_raster)

- **output_path**  
  Path where the final Output RTC raster path be saved.

- **output_db** *(default: `False`)*  
  If True, output in dB scale (10*log10). Default: False (linear)

- **block_size** *(default: `4096`)*  
  Size of processing blocks in pixels of raster for memory management.

- **verbose** *(default: `True`)*  
  If `True`, print verbose output during processing.

3- This function is responsible for getting the rotation matrix to convert the ECEF coordinate system to the ENU coordinate system at a given latitude and longitude. 

 ```python
ecef_to_enu_rotation_matrix(
    lat_deg: float,
    lon_deg: float,
    array_module)
 ```
#### Parameters

- **lat_deg**  
  Latitude in degrees 

- **lon_deg**  
  Longitude in degrees 

- **array_module**  *(default: None)*  
  The array module to use (e.g., `numpy` or `cupy`) , detects automatically if None.

4- This function applies the ENU to ECEF rotation matrix to a vector in ENU coordinates.

 ```python
ecef_to_enu_rotation_matrix(
    vector_ecef: float,
    scene_center_ecef: float,
    scene_geo ,
    is_position: bool = True)
 ```
#### Parameters

- **vector_ecef**  
  Vector in ECEF coordinates [X, Y, Z]

- **scene_center_ecef**  
  The origin of the ENU frame [X, Y, Z].

- **scene_geo**  
  Scene center geographic coordinates [Lat, Lon] in degrees

- **is_position** *(default: `True`)*  
  True if vector_ecef is a position (subtracts center), 
  False if it's a direction/velocity vector.

- **array_module**  *(default: None)*  
  The array module to use (e.g., `numpy` or `cupy`) , detects automatically if None.



5- This function converts ENU coordinates (meters) to Elevation Angle (degrees).

 ```python
ecef_to_enu_rotation_matrix(
    vector_ecef: float,
    scene_center_ecef: float,
    scene_geo ,
    is_position: bool = True)
 ```
#### Parameters

- **vector_ecef**  
  Vector in ECEF coordinates [X, Y, Z]

- **scene_center_ecef**  
  The origin of the ENU frame [X, Y, Z].

- **scene_geo**  
  Scene center geographic coordinates [Lat, Lon] in degrees

- **is_position** *(default: `True`)*  
  True if vector_ecef is a position (subtracts center), 
  False if it's a direction/velocity vector.

- **array_module**  *(default: None)*  
  The array module to use (e.g., `numpy` or `cupy`) , detects automatically if None.


6- This function converts ENU coordinates (meters) to Elevation Angle (degrees).

 ```python
enu_to_elevation_angle(
    east: float,
    north: float,
    up: float ,
    array_module=None)
 ```
#### Parameters

- **east** *(float)*  
  East component (meters)

- **north** *(float)*  
  North component (meters)

- **up** *(float)*  
  Up/Altitude component (meters)

- **array_module** *(default: None)*  
  The array module to use (e.g., `numpy` or `cupy`) , detects automatically if None.

- **array_module**  *(default: None)*  
  The array module to use (e.g., `numpy` or `cupy`) , detects automatically if None.


7- Converting a bounding box from a source CRS to a target CRS (default WGS84).

 ```python
transform_bbox(
    bbox: tuple[float, float, float, float],
    source_crs: str,
    target_crs: str = 'EPSG:4326'
)
 ```
#### Parameters

- **bbox** *(tuple[float, float, float, float])*  
  Bounding box in source CRS coordinates (min_x, min_y, max_x, max_y)

- **source_crs** *(str)*  
  Source coordinate reference system (e.g., 'EPSG:32633')

- **target_crs** *(str, default: 'EPSG:4326')*  
  Target coordinate reference system (e.g., 'EPSG:4326')


- **scene_center_ecef**  
  The origin of the ENU frame [X, Y, Z].

- **scene_geo**  
  Scene center geographic coordinates [Lat, Lon] in degrees

- **is_position** *(default: `True`)*  
  True if vector_ecef is a position (subtracts center), 
  False if it's a direction/velocity vector.

- **array_module**  *(default: None)*  
  The array module to use (e.g., `numpy` or `cupy`) , detects automatically if None.



### image_processing
- Purpose: TODO — describe filters, radiometric corrections, and preprocessing steps.
- Goals:
  - TODO
- Planned functions/classes:
  - TODO

### distortions
- Purpose: TODO — describe geometric distortion modelling and correction utilities.
- Goals:
  - TODO
- Planned functions/classes:
  - TODO

### io
- Purpose: TODO — describe input/output helpers, metadata readers/writers, and file converters.
- Goals:
  - TODO
- Planned functions/classes:
  - TODO

### ml
- Purpose: TODO — describe machine learning components: models, training utilities, inference helpers.
- Goals:
  - TODO
- Planned functions/classes:
  - TODO

### utils
- Purpose: TODO — general utilities and helpers (e.g., `dem_utils`).
- Goals:
  - TODO
- Planned functions/classes:
  - TODO

### Legacy
- Purpose: TODO — notes about legacy scripts and compatibility wrappers (see `Legacy/` folder).
- Migration plan:
  - TODO

### Docs
- Purpose: TODO — central documentation notes and additional specs (see `Docs/` folder).
  - TODO: link to Specs.md and other documentation pages

### Tests
- Purpose: TODO — testing plan and fixtures (see `tests/` folder).
  - TODO: describe test coverage goals

## Installation

Instructions to install the library (placeholder).

Recommended: create a virtual environment and install dependencies.

```bash
python -m venv .venv
source .venv/bin/activate  # or `.venv\\Scripts\\activate` on Windows
pip install -r requirements.txt
```

## Quick Start

Short examples and commands to get started. Add code snippets here as you implement modules.

## Usage Examples

- Example 1: TODO — simple DEM workflow
- Example 2: TODO — geometric transform + correction

## Contributing

Please add feature descriptions under each module and open issues/PRs for major changes. Include tests and update this README with usage examples.

## Roadmap

- Short-term: TODO
- Medium-term: TODO
- Long-term: TODO

## License

Add your license here (e.g., MIT). TODO: choose license.

---

Fill each module's section iteratively; keep this file as the central library documentation hub.
