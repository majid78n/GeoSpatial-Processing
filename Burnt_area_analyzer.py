"""
Burnt Area Quick Analyzer
A simplified Python library for quick burnt area analysis from satellite imagery using NBR.
"""

import rasterio
from rasterio.warp import transform_bounds, calculate_default_transform, reproject, Resampling
from rasterio.crs import CRS
import numpy as np
import matplotlib.pyplot as plt
import folium
from typing import Optional
from IPython.display import display
import os
import pyproj
os.environ["PROJ_LIB"] = pyproj.datadir.get_data_dir()

class QuickAnalyzer:
    """
    A simplified class for quick burnt area analysis.
    Features: image info, RGB visualization, NBR calculation, burnt area mapping, 
    area calculation, and location display on OpenStreetMap.
    """
    
    def __init__(self, filepath: str, nir_band: int = 4, swir_band: int = 5,
                 red_band: int = 1, green_band: int = 2, blue_band: int = 3):
        """
        Initialize the QuickAnalyzer.
        
        Args:
            filepath: Path to the satellite imagery file
            nir_band: Band number for Near-Infrared (default: 4)
            swir_band: Band number for Shortwave Infrared (default: 5)
            red_band: Band number for Red (default: 1)
            green_band: Band number for Green (default: 2)
            blue_band: Band number for Blue (default: 3)
        """
        self.filepath = filepath
        self.dataset = rasterio.open(filepath)
        self.nir_band = nir_band
        self.swir_band = swir_band
        self.red_band = red_band
        self.green_band = green_band
        self.blue_band = blue_band
    
    def show_info(self):
        """Display image information."""
        print("\n" + "=" * 60)
        print("IMAGE INFORMATION")
        print("=" * 60)
        print(f"File: {self.filepath}")
        print(f"Number of bands: {self.dataset.count}")
        print(f"Image size: {self.dataset.width} x {self.dataset.height} pixels")
        print(f"CRS: {self.dataset.crs}")
        print(f"Bounds: {self.dataset.bounds}")
        print(f"Resolution: {self.dataset.res}")
        if self.dataset.descriptions[0]:
            print(f"Band descriptions: {self.dataset.descriptions}")
        print("=" * 60 + "\n")
    
    @staticmethod
    def normalize_band(band: np.ndarray, percentile_range: tuple = (2, 98)) -> np.ndarray:
        """
        *** Percentiles make the normalization robust to outliers ***
        Normalize a band using percentile-based stretching.
        
        Args:
            band: Input band array
            percentile_range: Tuple of (lower, upper) percentile values
            
        Returns:
            Normalized band array with values between 0 and 1

        """
        band_valid = band[band > 0]
        
        if band_valid.size > 0:
            p_low, p_high = np.percentile(band_valid, percentile_range)
        else:
            p_low, p_high = band.min(), band.max()
        
        return np.interp(band, (p_low, p_high), (0, 1))
    
    def show_rgb_image(self, save_path: Optional[str] = None):
        """
        Display RGB composite of the satellite image.
        
        Args:
            save_path: Path to save the figure (optional)
        """
        print("Creating RGB composite...")
        
        red = self.dataset.read(self.red_band)
        green = self.dataset.read(self.green_band)
        blue = self.dataset.read(self.blue_band)
        
        rgb_composite = np.dstack((
            self.normalize_band(red),
            self.normalize_band(green),
            self.normalize_band(blue)
        ))
        
        plt.figure(figsize=(12, 10))
        plt.imshow(rgb_composite)
        plt.title("RGB Composite Image", fontsize=14, fontweight='bold')
        plt.axis('off')
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            print(f"RGB image saved to: {save_path}")
        
        plt.show()
        print("RGB visualization complete!")
    
    def show_location_map(self, save_path: Optional[str] = None):
        """
        Visualize the spatial footprint and center point of the raster dataset
        on an OpenStreetMap basemap (WGS84).
        """
        print("Creating location map...")

        # Get raster bounds and CRS
        bounds = self.dataset.bounds  # left, bottom, right, top in dataset CRS
        crs = self.dataset.crs        # coordinate reference system of the raster

        # Reproject bounds to WGS84 (EPSG:4326) if the raster is in a different CRS
        if crs and crs.to_string() != "EPSG:4326":
            bounds_wgs84 = transform_bounds(
                crs,                   # source CRS
                "EPSG:4326",           # target CRS
                bounds.left,
                bounds.bottom,
                bounds.right,
                bounds.top,
                densify_pts=21         # densify edges to approximate curved lines
            )
        else:
            bounds_wgs84 = (
                bounds.left,
                bounds.bottom,
                bounds.right,
                bounds.top
            )

        # Unpack bounds
        west, south, east, north = bounds_wgs84

        # Compute center coordinates of the raster footprint
        center_lat = (south + north) / 2
        center_lon = (west + east) / 2

        # Create a folium map centered at the raster center
        m = folium.Map(
            location=[center_lat, center_lon],
            zoom_start=10,       # initial zoom level
            tiles="OpenStreetMap" # basemap tiles
        )

        # Draw the raster bounding box as a rectangle
        folium.Rectangle(
            bounds=[[south, west], [north, east]],  # rectangle corners
            color="red",                           # border color
            fill=True,
            fill_opacity=0.2,                      # transparency of fill
            popup=f"Study Area<br>{self.filepath.split('/')[-1]}"  # popup text
        ).add_to(m)

        # Add a marker at the center of the raster
        folium.Marker(
            [center_lat, center_lon],
            popup=f"Center<br>Lat: {center_lat:.5f}<br>Lon: {center_lon:.5f}",
            icon=folium.Icon(color="red")
        ).add_to(m)
        
        # Optional: save the map to an HTML file if save_path is provided
        if save_path:
            # Ensure save_path includes a filename ending with .html
            if save_path.endswith("/") or save_path.endswith("\\"):
                save_path = os.path.join(save_path, "location_map.html")
            m.save(save_path)
            print(f"Location map saved to: {save_path}")

        print(f"Location (WGS84): Lat {center_lat:.5f}, Lon {center_lon:.5f}")
        print("Map created successfully.")

        return m  # Return the interactive map object for Jupyter display or further use
        
    def calculate_nbr(self) -> np.ndarray:
        """
        Calculate Normalized Burn Ratio (NBR).
            
        NBR = (NIR - SWIR) / (NIR + SWIR)
            
        Returns:
            NBR array
        """
        print("Calculating NBR...")
        nir = self.dataset.read(self.nir_band).astype(float)
        swir = self.dataset.read(self.swir_band).astype(float)
                
        with np.errstate(divide='ignore', invalid='ignore'):
            nbr = (nir - swir) / (nir + swir)
                
        print("NBR calculation complete!")
        return nbr
    
    def show_burnt_area_map(self, nbr: np.ndarray, threshold: float = -0.35, 
                           save_path: Optional[str] = None):
        """
        Visualize burnt areas on a map.
        
        Args:
            nbr: NBR array
            threshold: NBR threshold for burnt area detection
            save_path: Path to save the figure (optional)
        """
        print(f"Detecting burnt areas (NBR < {threshold})...")
        # Create a boolean mask of burnt areas: True where NBR is below the threshold
        burnt_mask = nbr < threshold

        # Create a figure with 2 side-by-side subplots
        fig, axes = plt.subplots(1, 2, figsize=(16, 7))

        # -------------------------------
        # Left subplot: display NBR values
        # -------------------------------
        # Use a diverging colormap to show variation from low to high NBR
        im1 = axes[0].imshow(nbr, cmap='RdYlGn', vmin=-1, vmax=1)
        axes[0].set_title('Normalized Burn Ratio (NBR)', fontsize=14, fontweight='bold')
        axes[0].axis('off')  # Hide axis ticks and labels
        # Add a colorbar for reference
        plt.colorbar(im1, ax=axes[0], fraction=0.046, pad=0.04, label='NBR Value')

        # -----------------------------------
        # Right subplot: visualize burnt areas
        # -----------------------------------
        # Create an RGBA image for visualization
        burnt_image = np.zeros((*burnt_mask.shape, 4))  # shape = (rows, cols, 4)

        # Set burnt pixels to red (R=1, G=0, B=0, alpha=0.8)
        burnt_image[burnt_mask] = [1, 0, 0, 0.8]

        # Set non-burnt pixels to light green with partial transparency
        burnt_image[~burnt_mask] = [0.2, 0.8, 0.2, 0.3]

        # Display the burnt area image
        axes[1].imshow(burnt_image)
        axes[1].set_title(f'Burnt Areas (NBR < {threshold})', fontsize=14, fontweight='bold')
        axes[1].axis('off')  # Hide axis ticks and labels

        # Adjust spacing between subplots
        plt.tight_layout()

        # Save the figure to a file if a save path is provided
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            print(f"Map saved to: {save_path}")

        # Show the figure interactively (useful in Jupyter)
        plt.show()
        print("Visualization complete!")
    
    def calculate_area(self, nbr: np.ndarray, threshold: float = -0.35) -> float:
        """
        Calculate total burnt area in km².
        
        Args:
            nbr: NBR array
            threshold: NBR threshold for burnt area detection
            
        Returns:
            Total burnt area in km²
        """
        print("Calculating burnt area...")
        burnt_mask = nbr < threshold
        
        # Count burnt pixels
        burnt_pixels = np.sum(burnt_mask)
        
        # Calculate pixel area
        pixel_width, pixel_height = self.dataset.res
        
        # Check if CRS is geographic (lat/lon in degrees)
        if self.dataset.crs and self.dataset.crs.is_geographic:
            # Convert degrees to meters
            bounds = self.dataset.bounds
            center_lat = (bounds.bottom + bounds.top) / 2
            
            # Conversion factors
            # Approximate conversion factors: meters per degree
            meters_per_deg_lat = 111320 # 1 degree latitude ≈ 111.32 km
            meters_per_deg_lon = 111320 * np.cos(np.radians(center_lat)) # Adjust for latitude
            
            # Convert pixel dimensions to meters
            pixel_width_m = abs(pixel_width * meters_per_deg_lon)
            pixel_height_m = abs(pixel_height * meters_per_deg_lat)
            
            # Area of one pixel in square meters
            pixel_area_m2 = pixel_width_m * pixel_height_m
            
            print(f"\nGeographic CRS detected (degrees)")
            print(f"Center latitude: {center_lat:.4f}°")
            print(f"Pixel size: {pixel_width:.6f}° x {pixel_height:.6f}°")
            print(f"Pixel size in meters: {pixel_width_m:.2f} x {pixel_height_m:.2f} m")
        else:
            # CRS is projected (already in meters)
            pixel_area_m2 = abs(pixel_width * pixel_height)
            print(f"\nProjected CRS detected")
            print(f"Pixel size: {pixel_width:.2f} x {pixel_height:.2f} m")
        
        # Calculate total area
        total_area_m2 = burnt_pixels * pixel_area_m2  # Total area in m²
        total_area_km2 = total_area_m2 / 1e6   # Convert to km²
        
        print(f"Burnt pixels: {burnt_pixels:,}")
        print(f"Pixel area: {pixel_area_m2:.2f} m²")
        print(f"Total burnt area: {total_area_km2:.2f} km²")
        
        return total_area_km2
    
    def analyze(self, threshold: float = -0.35, show_rgb: bool = True, 
                show_location: bool = True, save_outputs: bool = False,
                output_prefix: str = "output") -> float:
        """
        Run quick analysis with all features.
        
        Args:
            threshold: NBR threshold for burnt area detection
            show_rgb: Whether to display RGB composite
            show_location: Whether to display location map
            save_outputs: Whether to save all outputs
            output_prefix: Prefix for output files
            
        Returns:
            Total burnt area in km²
        """
        print("\n" + "="*60)
        print("QUICK BURNT AREA ANALYSIS")
        print("="*60)
        
        # Step 1: Show image info
        self.show_info()
        
        # Step 2: Show RGB composite (optional)
        if show_rgb:
            save_path = f"{output_prefix}_rgb.png" if save_outputs else None
            self.show_rgb_image(save_path=save_path)
            print()
        
        # Step 3: Show location map (optional)
        if show_location:
            save_path = f"{output_prefix}_location.html" if save_outputs else None
            m = self.show_location_map(save_path=save_path)#self.show_location_map(save_path=save_path)
            display(m)
            print()
        
        # Step 4: Calculate NBR
        nbr = self.calculate_nbr()
        
        # Step 5: Show burnt area map
        print()
        save_path = f"{output_prefix}_burnt_areas.png" if save_outputs else None
        self.show_burnt_area_map(nbr, threshold=threshold, save_path=save_path)
        
        # Step 6: Calculate and display area
        print()
        total_area = self.calculate_area(nbr, threshold=threshold)
        
        print("\n" + "="*60)
        print(f"ANALYSIS COMPLETE - Total Burnt Area: {total_area:.2f} km²")
        print("="*60 + "\n")
        
        return total_area
    
    def __del__(self):
        """Close the dataset when object is destroyed."""
        if hasattr(self, 'dataset'):
            self.dataset.close()


# Convenience function for quick analysis
def quick_analyze(filepath: str, threshold: float = -0.35, 
                  show_rgb: bool = True, show_location: bool = True,
                  save_outputs: bool = False, output_prefix: str = "output",
                  **kwargs) -> float:
    """
    Convenience function for quick burnt area analysis.
    
    Args:
        filepath: Path to satellite imagery file
        threshold: NBR threshold for burnt area detection
        show_rgb: Whether to display RGB composite
        show_location: Whether to display location map
        save_outputs: Whether to save all outputs
        output_prefix: Prefix for output files
        **kwargs: Additional arguments for QuickAnalyzer initialization (band numbers)
        
    Returns:
        Total burnt area in km²
    
    Example:
        area = quick_analyze('satellite_image.tif', threshold=-0.35, save_outputs=True)
    """
    analyzer = QuickAnalyzer(filepath, **kwargs)
    return analyzer.analyze(
        threshold=threshold,
        show_rgb=show_rgb,
        show_location=show_location,
        save_outputs=save_outputs,
        output_prefix=output_prefix
    )


if __name__ == "__main__":
    print("Burnt Area Quick Analyzer")
    print("\nQuick and simple burnt area analysis from satellite imagery")
    print("\nExample usage:")
    print("  analyzer = QuickAnalyzer('image.tif')")
    print("  area = analyzer.analyze()")
    print("\n  OR use the convenience function:")

    print("  area = quick_analyze('image.tif', save_outputs=True)")
