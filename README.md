# 🔥 Burnt Area Quick Analyzer

A streamlined Python library for rapid burnt area analysis from satellite imagery using the Normalized Burn Ratio (NBR) index.

[![Python 3.7+](https://img.shields.io/badge/python-3.7+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## 📋 Overview

Burnt Area Quick Analyzer is a simplified yet powerful tool designed for quick assessment of burnt areas from multispectral satellite imagery. Perfect for researchers, environmental scientists, and GIS analysts who need rapid fire damage assessments.

### Key Features

- 🖼️ **RGB Visualization**: Generate true-color composite images from multispectral data
- 📊 **NBR Calculation**: Automatic Normalized Burn Ratio computation
- 🗺️ **Interactive Maps**: Visualize burnt areas and display location on OpenStreetMap
- 📏 **Area Calculation**: Precise burnt area measurements in km²
- 🎯 **One-Line Analysis**: Run complete analysis with a single function call
- 💾 **Export Capabilities**: Save all outputs (images, maps, statistics)

## 🚀 Installation

### Prerequisites

```bash
conda install rasterio numpy matplotlib folium ipython
```

### Install from Source

```bash
git clone https://github.com/majid78n/GeoSpatial-Processing.git
cd burnt-area-analyzer
conda env create -f Requirements.yml
```

## 📖 Quick Start

### Basic Usage

```python
from burnt_area_analyzer import quick_analyze

# Analyze burnt area with one line
total_area = quick_analyze('satellite_image.tif', threshold=-0.35, save_outputs=True)
print(f"Total burnt area: {total_area:.2f} km²")
```
## Important Note on Band Selection

This library is designed primarily for **clipped Sentinel satellite images**.  
Band numbers for Sentinel are **different from Landsat, MODIS, or other satellite datasets**.

**To use this library with other images**, you must explicitly define the correct bands when creating a `QuickAnalyzer` object. For example:

```python
from burnt_area_analyzer import QuickAnalyzer

analyzer = QuickAnalyzer(
    "your_image.tif",
    nir_band=4,    # Near-Infrared band
    swir_band=5,   # Shortwave Infrared band
    red_band=1,
    green_band=2,
    blue_band=3
)
```
⚠️ Warning: Incorrect band assignment can lead to wrong RGB visualization or invalid NBR/burnt area calculations. Always check your image’s band order before analysis.

### Advanced Usage with Class

```python
from burnt_area_analyzer import QuickAnalyzer

# Initialize analyzer
analyzer = QuickAnalyzer(
    filepath='satellite_image.tif',
    nir_band=4,    # Near-Infrared band
    swir_band=5,   # Shortwave Infrared band
    red_band=1,
    green_band=2,
    blue_band=3
)

# Show image information
analyzer.show_info()

# Display RGB composite
analyzer.show_rgb_image(save_path='rgb_output.png')

# Show location on map
map_object = analyzer.show_location_map(save_path='location.html')

# Calculate NBR
nbr = analyzer.calculate_nbr()

# Visualize burnt areas
analyzer.show_burnt_area_map(nbr, threshold=-0.35, save_path='burnt_areas.png')

# Calculate total burnt area
total_area = analyzer.calculate_area(nbr, threshold=-0.35)
print(f"Burnt area: {total_area:.2f} km²")
```

### Complete Analysis Pipeline

```python
from burnt_area_analyzer import QuickAnalyzer

analyzer = QuickAnalyzer('satellite_image.tif')

# Run complete analysis with all features
area = analyzer.analyze(
    threshold=-0.35,        # NBR threshold for burnt detection
    show_rgb=True,          # Display RGB composite
    show_location=True,     # Show location map
    save_outputs=True,      # Save all outputs
    output_prefix='fire_2024'  # Prefix for output files
)
```

## 🎯 How It Works

### Normalized Burn Ratio (NBR)

The library uses the NBR index to identify burnt areas:

```
NBR = (NIR - SWIR) / (NIR + SWIR)
```

Where:
- **NIR**: Near-Infrared band (typically band 4 in Landsat)
- **SWIR**: Shortwave Infrared band (typically band 5 in Landsat)

**NBR Values Interpretation:**
- **High NBR** (0.1 to 1.0): Healthy vegetation
- **Low NBR** (-0.1 to 0.1): Bare soil or recent burns
- **Very Low NBR** (-1.0 to -0.1): Severe burns

Default threshold: **-0.35** (adjustable based on your needs)

## 📊 Outputs

### 1. Image Information
```
==============================================================
IMAGE INFORMATION
==============================================================
File: satellite_image.tif
Number of bands: 7
Image size: 7801 x 7651 pixels
CRS: EPSG:32633
Bounds: BoundingBox(left=300000, bottom=4500000, right=540000, top=4730000)
Resolution: (30.0, 30.0)
==============================================================
```

### 2. RGB Composite Image
- True-color visualization of the study area
- Percentile-based normalization for optimal contrast

### 3. Location Map
- Interactive OpenStreetMap with study area boundary
- Center point marker with coordinates
- HTML export for sharing

### 4. Burnt Area Visualization
- Side-by-side comparison: NBR values vs. Burnt areas
- Color-coded burnt regions (red) and healthy areas (green)

### 5. Area Statistics
```
Burnt pixels: 125,432
Pixel area: 900.00 m²
Total burnt area: 112.89 km²
```

## 🛠️ API Reference

### `QuickAnalyzer` Class

#### Constructor
```python
QuickAnalyzer(filepath, nir_band=4, swir_band=5, red_band=1, green_band=2, blue_band=3)
```

#### Methods

| Method | Description | Returns |
|--------|-------------|---------|
| `show_info()` | Display image metadata | None |
| `show_rgb_image(save_path=None)` | Display RGB composite | None |
| `show_location_map(save_path=None)` | Show location on map | folium.Map |
| `calculate_nbr()` | Calculate NBR index | np.ndarray |
| `show_burnt_area_map(nbr, threshold, save_path=None)` | Visualize burnt areas | None |
| `calculate_area(nbr, threshold)` | Calculate burnt area | float (km²) |
| `analyze(threshold, show_rgb, show_location, save_outputs, output_prefix)` | Run complete analysis | float (km²) |

### `quick_analyze()` Function

```python
quick_analyze(
    filepath,
    threshold=-0.35,
    show_rgb=True,
    show_location=True,
    save_outputs=False,
    output_prefix='output',
    **kwargs
)
```

**Parameters:**
- `filepath` (str): Path to satellite imagery file
- `threshold` (float): NBR threshold for burnt detection
- `show_rgb` (bool): Display RGB composite
- `show_location` (bool): Show location map
- `save_outputs` (bool): Save all outputs to files
- `output_prefix` (str): Prefix for output filenames
- `**kwargs`: Additional band configuration (nir_band, swir_band, etc.)

## 📁 Supported Data Formats

The library works with any raster format supported by `rasterio`:

- ✅ GeoTIFF (.tif, .tiff)
- ✅ HDF (.hdf, .h5)
- ✅ NetCDF (.nc)
- ✅ JP2 (.jp2)
- ✅ IMG (.img)
- ✅ And more...

### Satellite Imagery Compatibility

- **Landsat 8/9**: Bands 5 (NIR), 7 (SWIR2)
- **Sentinel-2**: Bands 8 (NIR), 12 (SWIR2)
- **Landsat 4/5/7**: Bands 4 (NIR), 7 (SWIR2)

## 🌍 Coordinate Reference Systems

The library automatically handles:
- Geographic CRS (lat/lon in degrees) - converts to meters for area calculation
- Projected CRS (already in meters) - uses directly
- Automatic reprojection to WGS84 for location maps

## 💡 Examples

### Example 1: Landsat 8 Analysis

```python
from burnt_area_analyzer import QuickAnalyzer

# Landsat 8 band configuration
analyzer = QuickAnalyzer(
    'LC08_scene.tif',
    nir_band=5,   # Landsat 8 NIR
    swir_band=7,  # Landsat 8 SWIR2
    red_band=4,
    green_band=3,
    blue_band=2
)

area = analyzer.analyze(threshold=-0.3)
```

### Example 2: Sentinel-2 Analysis

```python
from burnt_area_analyzer import QuickAnalyzer

# Sentinel-2 band configuration
analyzer = QuickAnalyzer(
    'S2A_scene.tif',
    nir_band=8,   # Sentinel-2 NIR
    swir_band=12, # Sentinel-2 SWIR2
    red_band=4,
    green_band=3,
    blue_band=2
)

area = analyzer.analyze(threshold=-0.35, save_outputs=True)
```

### Example 3: Custom Threshold Sensitivity Test

```python
from burnt_area_analyzer import QuickAnalyzer

analyzer = QuickAnalyzer('wildfire_scene.tif')
nbr = analyzer.calculate_nbr()

# Test different thresholds
thresholds = [-0.2, -0.3, -0.4, -0.5]

for threshold in thresholds:
    area = analyzer.calculate_area(nbr, threshold=threshold)
    print(f"Threshold {threshold}: {area:.2f} km²")
```

## 🔧 Troubleshooting

### Common Issues

**Issue: "Band index out of range"**
```python
# Solution: Check your band configuration
analyzer = QuickAnalyzer('image.tif', nir_band=5, swir_band=7)  # Adjust band numbers
```

**Issue: "Cannot display map in terminal"**
```python
# Solution: Use Jupyter Notebook or save to HTML
m = analyzer.show_location_map(save_path='map.html')
# Open map.html in your browser
```

**Issue: "Area calculation seems incorrect"**
```python
# Solution: Verify CRS and check pixel resolution
analyzer.show_info()  # Check resolution and CRS
```

## 📚 Scientific Background

### References

1. **Key, C. H., & Benson, N. C. (2006)**. "Landscape Assessment: Ground measure of severity, the Composite Burn Index; and Remote sensing of severity, the Normalized Burn Ratio." FIREMON: Fire Effects Monitoring and Inventory System.

2. **García, M. J. L., & Caselles, V. (1991)**. "Mapping burns and natural reforestation using Thematic Mapper data." Geocarto International, 6(1), 31-37.

### NBR Burn Severity Classification

| NBR Range | Burn Severity |
|-----------|---------------|
| > 0.1 | Unburned |
| 0.1 to -0.1 | Low severity |
| -0.1 to -0.27 | Moderate-low |
| -0.27 to -0.44 | Moderate-high |
| < -0.44 | High severity |


## 📝 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 👥 Authors

- Rohollah Naeijian - *M.Sc. student in Geoinformatics engineering at PoliMi*

## 📧 Contact

- GitHub: [@majid78n](https://github.com/majid78n)
- Email: rohollah.naeijian@mail.polimi.it
<img src="https://upload.wikimedia.org/wikipedia/it/b/be/Logo_Politecnico_Milano.png" alt="PoliMi Logo" width="100"/> 


