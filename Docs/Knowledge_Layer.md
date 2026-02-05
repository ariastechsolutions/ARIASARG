# SARPRO Knowledge Layer
## Mathematical Formulas and Equations Reference

**Version:** 1.0
**Last Updated:** February 5, 2026
**Purpose:** Comprehensive documentation of all mathematical formulas, algorithms, and computational methods used in SARPRO library.

---

## Table of Contents

1. [RPC (Rational Polynomial Coefficients)](#1-rpc-rational-polynomial-coefficients)
2. [Coordinate Transformations](#2-coordinate-transformations)
3. [SAR Radiometric Corrections](#3-sar-radiometric-corrections)
4. [Lee Speckle Filter](#4-lee-speckle-filter)
5. [DEM Operations](#5-dem-operations)
6. [Inverse RPC and Newton-Raphson](#6-inverse-rpc-and-newton-raphson)

---

## 1. RPC (Rational Polynomial Coefficients)

### 1.1 Coordinate Normalization

**Purpose:** Normalize geographic and image coordinates to improve numerical stability in RPC calculations.

#### Geographic Coordinate Normalization

```
L = (lon - lon_off) / lon_scale
P = (lat - lat_off) / lat_scale
H = (height - height_off) / height_scale
```

Where:
- `L` = Normalized longitude
- `P` = Normalized latitude
- `H` = Normalized height
- `lon_off`, `lat_off`, `height_off` = Offset values from RPC metadata
- `lon_scale`, `lat_scale`, `height_scale` = Scale values from RPC metadata

**Location:** `geometry/rpc_utils.py::normalize_coordinates()`

---

#### Image Coordinate Denormalization

```
row = row_n × line_scale + line_off
col = col_n × samp_scale + samp_off
```

Where:
- `row_n`, `col_n` = Normalized image coordinates
- `row`, `col` = Pixel coordinates
- `line_off`, `samp_off` = Image offset values
- `line_scale`, `samp_scale` = Image scale values

**Location:** `geometry/rpc_utils.py::denormalize_image_coordinates()`

---

### 1.2 RPC Polynomial Evaluation

**Purpose:** Evaluate the 20-term RPC polynomial used in satellite image georeferencing.

#### Standard 20-Term RPC Polynomial

```
f(P, L, H) = c₀ + c₁L + c₂P + c₃H + c₄LP + c₅LH + c₆PH + c₇L² + c₈P² + c₉H²
           + c₁₀PLH + c₁₁L³ + c₁₂LP² + c₁₃LH² + c₁₄L²P + c₁₅P³ + c₁₆PH²
           + c₁₇L²H + c₁₈P²H + c₁₉H³
```

Where:
- `P` = Normalized latitude
- `L` = Normalized longitude
- `H` = Normalized height
- `c₀...c₁₉` = RPC polynomial coefficients

**Terms:**
1. Constant: `1`
2. Linear: `L`, `P`, `H`
3. Cross products: `LP`, `LH`, `PH`
4. Quadratic: `L²`, `P²`, `H²`
5. Cubic cross: `PLH`, `L³`, `LP²`, `LH²`, `L²P`, `P³`, `PH²`, `L²H`, `P²H`, `H³`

**Location:** `geometry/rpc_polynomial.py::evaluate_rpc_polynomial()`

---

#### RPC Rational Function

**Forward RPC Transformation (Ground → Image):**

```
row_normalized = P_num(P, L, H) / P_den(P, L, H)
col_normalized = S_num(P, L, H) / S_den(P, L, H)
```

Where:
- `P_num`, `P_den` = Line (row) numerator and denominator polynomials
- `S_num`, `S_den` = Sample (column) numerator and denominator polynomials

Each polynomial follows the 20-term standard form.

**Location:** `geometry/rpc_polynomial.py::evaluate_rpc_rational_function()`

---

## 2. Coordinate Transformations

### 2.1 ECEF to ENU Rotation Matrix

**Purpose:** Transform Earth-Centered, Earth-Fixed (ECEF) coordinates to local East-North-Up (ENU) coordinates.

#### Rotation Matrix

```
R = [ -sin(λ)              cos(λ)               0        ]
    [ -sin(φ)cos(λ)       -sin(φ)sin(λ)       cos(φ)    ]
    [  cos(φ)cos(λ)        cos(φ)sin(λ)       sin(φ)    ]
```

Where:
- `φ` = Latitude (radians)
- `λ` = Longitude (radians)
- `R` = 3×3 rotation matrix

**Transformation:**
```
[E]       [X_ecef - X_center]
[N] = R · [Y_ecef - Y_center]
[U]       [Z_ecef - Z_center]
```

Where:
- `E`, `N`, `U` = East, North, Up components
- `X_ecef`, `Y_ecef`, `Z_ecef` = ECEF coordinates
- `X_center`, `Y_center`, `Z_center` = Local origin in ECEF

**Location:** `geometry/transforms.py::ecef_to_enu_rotation_matrix()`

---

### 2.2 ENU to Elevation Angle

**Purpose:** Convert local ENU coordinates to elevation angle (altitude angle from horizon).

```
ground_dist = √(E² + N²)
elevation_rad = arctan2(U, ground_dist)
elevation_deg = elevation_rad × (180/π)
```

Where:
- `E`, `N`, `U` = East, North, Up components (meters)
- `ground_dist` = Horizontal distance
- `elevation_deg` = Elevation angle in degrees (0° = horizon, 90° = zenith)

**Location:** `geometry/transforms.py::enu_to_elevation_angle()`

---

### 2.3 Affine Transformation

**Purpose:** Transform map coordinates to pixel coordinates using affine transformation.

#### Forward Affine Transform

```
x_map = a × col + b × row + c
y_map = d × col + e × row + f
```

#### Inverse Affine Transform (Map → Pixel)

```
col = a⁻¹ × x_map + b⁻¹ × y_map + c⁻¹
row = d⁻¹ × x_map + e⁻¹ × y_map + f⁻¹
```

Where:
- `(x_map, y_map)` = Map coordinates (meters or degrees)
- `(col, row)` = Pixel coordinates
- `a, b, c, d, e, f` = Affine transformation parameters
- `a⁻¹, b⁻¹, c⁻¹, d⁻¹, e⁻¹, f⁻¹` = Inverse affine parameters

**Location:** `geometry/transforms.py::apply_affine_transform()`

---

## 3. SAR Radiometric Corrections

### 3.1 RTC Geometric Factor

**Purpose:** Compute the geometric normalization factor for Radiometric Terrain Correction (RTC) of SAR imagery.

```
γ = sin(θ_local) / sin(θ_ellipsoid)
```

Where:
- `γ` = Geometric normalization factor
- `θ_local` = Local incidence angle (radians) - angle between radar beam and terrain normal
- `θ_ellipsoid` = Ellipsoid incidence angle (degrees) - angle between radar beam and ellipsoid normal

**Implementation details:**
- `θ_ellipsoid` is converted from degrees to radians
- Division by zero is prevented by clamping `sin(θ_ellipsoid)` to minimum 1e-6
- Negative values are clamped to 0.0

**Location:** `geometry/transforms.py::calculate_geometric_factor()`

**Reference:**
Small, D. (2011). *Flattening Gamma: Radiometric Terrain Correction for SAR Imagery*.
IEEE Transactions on Geoscience and Remote Sensing, 49(8), 3081-3093.

---

### 3.2 RTC Correction Application

**Purpose:** Apply radiometric terrain correction to SAR digital numbers.

#### Linear Power Scale

```
σ_raw = DN² / K
σ_rtc = σ_raw × γ
```

#### Logarithmic (dB) Scale

```
σ_rtc_dB = 10 × log₁₀(σ_rtc)
```

Where:
- `DN` = Digital Number (raw SAR pixel value)
- `K` = Sensor-specific calibration constant
- `σ_raw` = Raw backscatter coefficient
- `σ_rtc` = RTC-corrected backscatter coefficient
- `γ` = Geometric factor (from equation 3.1)
- `σ_rtc_dB` = RTC-corrected backscatter in decibels

**Special Cases:**
- Negative values are clamped to 0.0 before dB conversion
- For dB conversion, values ≤ 0 are set to -99.0 dB

**Location:** `geometry/transforms.py::apply_rtc_correction()`

---

## 4. Lee Speckle Filter

**Purpose:** Adaptive filter to reduce speckle noise in SAR imagery while preserving edges and features.

### 4.1 Noise Coefficient of Variation

```
Cu² = 1 / ENL
```

Where:
- `Cu` = Coefficient of variation of noise
- `ENL` = Equivalent Number of Looks (1.0 for single-look, higher for multi-looked)

---

### 4.2 Local Statistics

#### Local Mean

```
μ_local = (1/N) Σ I(i,j)  for all (i,j) in window
```

#### Local Variance

```
Var_local = E[I²] - (E[I])²
Var_local = μ_sqr - μ_local²
```

Where:
- `μ_local` = Local mean intensity
- `μ_sqr` = Mean of squared intensities
- `Var_local` = Local variance
- `I(i,j)` = Intensity at pixel (i,j)
- `N` = Number of pixels in window

---

### 4.3 Lee Weighting Factor

```
speckle_var = μ_local² × Cu²
K = 1 - (speckle_var / Var_local)
K_clipped = clip(K, 0, 1)
```

Where:
- `speckle_var` = Theoretical variance due to speckle
- `K` = Lee weighting factor
- `K_clipped` = Weight clipped to [0, 1] range

**Interpretation:**
- `K ≈ 0` in homogeneous areas → strong smoothing
- `K ≈ 1` near edges → preserve original values

---

### 4.4 Lee Filter Output

```
I_filtered = μ_local + K × (I_original - μ_local)
```

Where:
- `I_filtered` = Filtered pixel value
- `I_original` = Original pixel value
- `μ_local` = Local mean
- `K` = Lee weighting factor

**Location:** `image_processing/filters.py::lee_filter_gpu()`

**Reference:**
Lee, J.-S. (1980). *Digital image enhancement and noise filtering by use of local statistics*.
IEEE Transactions on Pattern Analysis and Machine Intelligence, PAMI-2(2), 165-168.

---

## 5. DEM Operations

### 5.1 DEM Sampling with Affine Transform

**Purpose:** Sample DEM elevation at arbitrary geographic coordinates using inverse affine transformation.

#### Geographic to Pixel Coordinates

```
dem_col = a × lon + b × lat + c
dem_row = d × lon + e × lat + f
```

Where:
- `lon`, `lat` = Geographic coordinates (degrees)
- `dem_col`, `dem_row` = DEM pixel coordinates
- `a, b, c, d, e, f` = Inverse affine transformation coefficients

#### Bilinear Interpolation

The elevation is sampled using bilinear interpolation (order=1) via `map_coordinates`:

```
elevation = bilinear_interpolate(DEM, [dem_row, dem_col])
```

**Location:** `dem/sampling.py::sample_dem_at_coordinates()`

---

## 6. Inverse RPC and Newton-Raphson

**Purpose:** Solve the inverse RPC problem (image coordinates → ground coordinates) using iterative Newton-Raphson method.

### 6.1 Problem Formulation

Given target image coordinates `(row_target, col_target)`, find ground coordinates `(lon, lat)` such that:

```
forward_RPC(lon, lat, h) = (row_target, col_target)
```

This is a non-linear system with no closed-form solution, requiring iterative methods.

---

### 6.2 Newton-Raphson Iteration

#### Residual Calculation

```
r_err = row_target - row_estimated
c_err = col_target - col_estimated
```

Where:
- `row_estimated`, `col_estimated` = Forward RPC projection of current estimate

---

#### Jacobian Matrix (Finite Differences)

```
∂r/∂P ≈ (r(P + δ) - r(P)) / δ
∂c/∂P ≈ (c(P + δ) - c(P)) / δ
∂r/∂L ≈ (r(L + δ) - r(L)) / δ
∂c/∂L ≈ (c(L + δ) - c(L)) / δ
```

Where:
- `P` = Normalized latitude
- `L` = Normalized longitude
- `δ` = Small perturbation (typically 1e-4)
- `r`, `c` = Row and column from forward RPC

**Jacobian Matrix:**
```
J = [ ∂r/∂P   ∂r/∂L ]
    [ ∂c/∂P   ∂c/∂L ]
```

**Location:** `geometry/rpc_jacobian.py::compute_jacobian_finite_difference()`

---

#### Linear System Solution (Cramer's Rule)

Solve: `J · [dP, dL]ᵀ = [r_err, c_err]ᵀ`

```
det = (∂r/∂P × ∂c/∂L) - (∂r/∂L × ∂c/∂P)

dP = (∂c/∂L × r_err - ∂r/∂L × c_err) / det
dL = (∂r/∂P × c_err - ∂c/∂P × r_err) / det
```

**Update step:**
```
P_new = P_old + dP
L_new = L_old + dL
```

**Location:** `geometry/rpc_jacobian.py::solve_newton_step()`

---

#### Convergence

The iteration continues for a fixed number of steps (typically 5-10) or until:

```
||[dP, dL]|| < tolerance
```

**Division by zero prevention:**
- If `|det| < 1e-10`, set `det = 1e-10`

---

### 6.3 DEM Integration

When DEM is available, height is updated at each iteration:

```
1. Denormalize current estimate: (lon, lat) = denormalize(P, L)
2. Sample DEM: h = sample_DEM(lon, lat)
3. Use updated height in forward RPC for next iteration
```

This refines the solution by accounting for actual terrain elevation.

---

## Appendix A: Notation and Conventions

### Coordinate Systems

| Symbol | Description | Units |
|--------|-------------|-------|
| `(lon, lat, h)` | Geographic coordinates | degrees, degrees, meters |
| `(P, L, H)` | Normalized geographic coordinates | dimensionless |
| `(row, col)` | Image pixel coordinates | pixels |
| `(X, Y, Z)` | ECEF coordinates | meters |
| `(E, N, U)` | ENU local coordinates | meters |

### SAR Parameters

| Symbol | Description | Units |
|--------|-------------|-------|
| `DN` | Digital Number (raw pixel value) | dimensionless |
| `σ` | Backscatter coefficient | linear or dB |
| `K` | Calibration constant | sensor-specific |
| `γ` | Geometric factor | dimensionless |
| `θ_local` | Local incidence angle | radians |
| `θ_ellipsoid` | Ellipsoid incidence angle | degrees |
| `ENL` | Equivalent Number of Looks | dimensionless |

---

## Appendix B: Algorithm Complexity

| Algorithm | Time Complexity | Space Complexity | Notes |
|-----------|----------------|------------------|-------|
| RPC Forward | O(1) per point | O(1) | 20-term polynomial evaluation |
| RPC Inverse | O(n × k) | O(1) | n iterations, k = RPC forward cost |
| Lee Filter | O(w² × N) | O(N) | w = window size, N = image pixels |
| DEM Sampling | O(1) per point | O(1) | Bilinear interpolation |
| ECEF to ENU | O(1) per point | O(1) | Matrix multiplication |

---

## Appendix C: Numerical Stability

### Safeguards Implemented

1. **Division by Zero Prevention:**
   - Clamping denominators to minimum threshold (typically 1e-6 or 1e-10)
   - Used in: RTC geometric factor, Newton-Raphson det

2. **Range Clipping:**
   - Lee filter weights clipped to [0, 1]
   - RTC geometric factor clipped to [0, ∞)

3. **Logarithm Domain Protection:**
   - Setting `log(x ≤ 0) = -99.0` for dB conversion

4. **Coordinate Normalization:**
   - RPC uses normalized coordinates to prevent overflow/underflow

---

## Version History

| Version | Date | Changes |
|---------|------|---------|
| 1.0 | 2026-02-05 | Initial knowledge layer documentation |

---

## References

1. Small, D. (2011). *Flattening Gamma: Radiometric Terrain Correction for SAR Imagery*. IEEE Transactions on Geoscience and Remote Sensing, 49(8), 3081-3093.

2. Lee, J.-S. (1980). *Digital image enhancement and noise filtering by use of local statistics*. IEEE Transactions on Pattern Analysis and Machine Intelligence, PAMI-2(2), 165-168.

3. NGA (2014). *NITF RPC00B specification*. National Geospatial-Intelligence Agency Standard NGA.STND.0014.

4. Farrell, J. A. (2008). *Aided Navigation: GPS with High Rate Sensors*. McGraw-Hill Professional.

---

**End of Knowledge Layer Documentation**
