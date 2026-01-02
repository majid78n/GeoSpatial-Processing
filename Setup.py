"""
Setup script for Burnt Area Analyzer library
"""

from setuptools import setup, find_packages
import os

# Read README for long description (optional - skip if README doesn't exist yet)
long_description = "A Python library for detecting and analyzing burnt areas from satellite imagery using NBR"
if os.path.exists("README.md"):
    with open("README.md", "r", encoding="utf-8") as fh:
        long_description = fh.read()

setup(
    name="burnt-area-analyzer",
    version="1.0.0",
    author="Majid",  # Replace with your name
    author_email="Rohollah.naeijian@mail.polimi.it",  # Replace with your email
    description="A Python library for detecting and analyzing burnt areas from satellite imagery",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/yourusername/burnt-area-analyzer",  # Replace with your GitHub URL
    py_modules=["burnt_area_analyzer"],  # This tells setup.py about your main file
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: GIS",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
    ],
    python_requires=">=3.8",
    install_requires=[
        "rasterio>=1.3.0",
        "numpy>=1.21.0",
        "matplotlib>=3.5.0",
        "geopandas>=0.12.0",
        "shapely>=2.0.0",
        "pandas>=1.3.0",
    ],
)