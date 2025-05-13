# Green-Space-Accessibility

This is a Python 3.9+ package for measuring accessibility to green space entrances, designed for urban planners and researchers analysing geospatial data for applications such as hedonic price models. The package uses OpenStreetMap (OSM) data via `osmnx` to identify walking networks and green space boundaries, calculating the total area of accessible green spaces within a specified walking time (default: 15 minutes) from the centroid of input zones (e.g., Traffic Analysis Zones, TAZs).

* The current version calculates accessibility by:
  * Downloading OSM walking networks and green space boundaries (parks, recreation grounds, woods, scrub, grass) within a buffered area around input zones.
  * Identifying green space entrances as intersections between walking networks and green space boundaries.
  * Creating isochrone polygons (catchments) based on walking time from zone centroids.
  * Summing the area of unique green spaces accessible within each catchment.
* Unlike traditional methods that use green space centroids, this package measures accessibility to entrances, providing a more accurate representation of how people access parks.
* The tool is flexible, accepting input geometries such as TAZs or other polygons with a unique identifier.

### Planned Features
* Support for additional accessibility metrics, such as the number of green spaces or distance to nearest entrance.
* Integration with hedonic price model workflows for direct application in property value analysis.
* A QGIS plugin to make the tool accessible to GIS users.

These functionalities will be available as a Python library, with potential extensions for GIS software integration.

For reference, this work builds on research conducted with Peter Nunns, as documented in:
Nunns, P., Allpress, J., & Balderston, K. (2016). *How do Aucklanders value their parks? A hedonic analysis of the impact of proximity to open space on residential property values.* Auckland Council Technical Report, TR2016/031. Available at: [https://knowledgeauckland.org.nz/media/1292/tr2016-031-how-do-aucklanders-value-their-parks.pdf](https://knowledgeauckland.org.nz/media/1292/tr2016-031-how-do-aucklanders-value-their-parks.pdf)

## Installation

Install with pip via:
```bash
pip install git+https://github.com/saeidadli/green-space-accessibility
```

Dependencies:
* `geopandas`
* `networkx`
* `osmnx`
* `pandas`
* `shapely`

Install dependencies manually with:
```bash
pip install geopandas networkx osmnx pandas shapely
```

Note: Ensure you have a compatible version of `gdal` installed for `geopandas`. On some systems, you may need to install it separately (e.g., `libgdal-dev` on Ubuntu).

## Usage

See the example Jupyter notebook at `examples/access2green_example.ipynb` for a step-by-step guide on:
* Loading TAZ data (e.g., from a GeoJSON file).
* Running the `access2green` function to calculate green space accessibility.
* Saving the output to a shapefile.

Basic example:
```python
import geopandas as gpd
import access2green as ag

# Load TAZ data
tazs = gpd.read_file("data/auckland.geojson")

# Calculate accessibility to green spaces
gdf = ag.access2green(
    tazs=tazs,
    uid="SA12020_V1",  # Column with unique ID for each zone
    buffer_distance_meters=1000,  # Buffer distance for OSM download
    trip_time=15,  # Walking time in minutes
    travel_speed=4.5  # Walking speed in km/h
)

# Save output
gdf.to_file("data/access2green.shp")
```

The output GeoDataFrame includes the original TAZ geometries with an additional `total_green_area` column, representing the sum of accessible green space areas (in square meters) within the walking catchment.

## Notes

* This package is experimental and relies on OpenStreetMap data, which may vary in quality depending on the region.
* The input GeoDataFrame must have a WGS 84 CRS (EPSG:4326). The function automatically projects to a suitable UTM CRS for calculations.
* No automated tests are implemented yet.
* Contributions, feedback, and suggestions are warmly welcomed.

## Authors

* Saeid Adli, 2025/05/13