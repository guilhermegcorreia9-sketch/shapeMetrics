<!-- README.md -->

# shapeMetrics

<!-- badges: start -->

[![R-CMD-check](https://github.com/guilhermegcorreia9-sketch/shapeMetrics/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/guilhermegcorreia9-sketch/shapeMetrics/actions/workflows/R-CMD-check.yaml) [![CRAN status](https://www.r-pkg.org/badges/version/shapeMetrics)](https://CRAN.R-project.org/package=shapeMetrics) [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

<!-- badges: end -->

<img src="man/figures/shapeMetrics_logo.png" align="right" height="139"/>

**shapeMetrics** is an R package that provides a comprehensive set of **shape and compactness indices** for polygon geometries (`sf` objects). It implements over 30 metrics inspired by the `vectormetrics` and `redistmetrics` packages, with built-in **parallel processing** for large datasets.

------------------------------------------------------------------------

## Installation

### From GitHub (development version)

``` r
# Install devtools if not already installed
install.packages("devtools")

# Install shapeMetrics from GitHub
devtools::install_github("guilhermegcorreia9-sketch/shapeMetrics")
```

## Quick Example

``` r
library(sf)
library(shapeMetrics)

# Load example polygons (North Carolina counties)
nc <- st_read(system.file("shape/nc.shp", package = "sf"))

# Calculate a single metric (Polsby-Popper) with 2 cores
nc <- polsby_popper(nc, ncores = 2)

# Calculate all metrics at once
nc <- calc_multiple_metrics(nc, ncores = 2)

# View results
head(nc[, c("NAME", "polsby_popper", "shape_index", "compactness")])
```

## Features

- **30+ shape and compactness metrics** – including Polsby-Popper, Reock, Schwartzberg, Boyce-Clark, and many more.
- **Parallel processing** – each metric is computed in parallel across polygons using `future.apply`.
- **Automatic handling of geographic CRS** – disables spherical geometry (s2) automatically for metrics that require planar calculations, and restores it afterward.
- **Robust error handling** – emits informative warnings when metrics fail (e.g., invalid geometry, sampling issues).
- **Works with `sf` objects** – all functions accept and return `sf` objects with new columns.

## List of Implemented Metrics

### Shape Metrics (vectormetrics style)

| Metric | Function | Description |
|----|----|----|
| Circularity | `circularity()` | A/(P²/4π) |
| Compactness (Form Factor) | `compactness()` | 4πA/P² |
| Convexity | `convexity()` | P_hull/P |
| Elongation | `elongation()` | Minor axis / Major axis |
| Fractal Dimension | `fractal_dimension()` | 2·log(P)/log(A) |
| Rectangularity (MABR) | `rectangularity()` | A/A_MABR |
| Roughness | `roughness()` | Shape Index − 1 |
| Shape Index | `shape_index()` | P/(2√(πA)) |
| Solidity (Convex Hull Compactness) | `convex_hull_compact()` | A/A_convex hull |
| Sphericity | `sphericity()` | r_inscribed/r_circumscribed |

### Compactness Metrics (redistmetrics style)

| Metric        | Function          | Description                       |
|---------------|-------------------|-----------------------------------|
| Polsby-Popper | `polsby_popper()` | 4πA/P²                            |
| Schwartzberg  | `schwartzberg()`  | 1/(P/(2π√(A/π)))                  |
| Reock         | `reock()`         | A/A_circumscribing circle         |
| Box-Reock     | `box_reock()`     | A/A_bounding rectangle            |
| Boyce-Clark   | `boyce_clark()`   | Radial distances from centroid    |
| Skew          | `skew()`          | A_inscribed/A_circumscribed       |
| Girth         | `girth()`         | r_inscribed/r_equal-area          |
| Detour        | `detour()`        | P_equal-area circle/P_convex hull |

### Additional Indices

| Metric | Function | Description |
|----|----|----|
| Cohesion | `cohesion()` | Mean squared distance ratio |
| Exchange | `exchange()` | Proportion of area inside equal-area circle |
| Proximity | `proximity()` | Mean distance ratio |
| Range Index | `range_index()` | Diameter ratio |
| Equivalent Rectangular | `equiv_rectangular()` | P_square/P |
| Squareness | `squareness()` | Same as equivalent rectangular |
| Perimeter Index | `perimeter_index()` | 1/shape_index |
| Perimeter-Area Ratio | `perimeter_area_ratio()` | P/A |

### Geometric Attributes

| Function                  | Description                        |
|---------------------------|------------------------------------|
| `area()`                  | Area of polygons                   |
| `perimeter()`             | Perimeter of polygons              |
| `convex_hull_area()`      | Area of convex hull                |
| `convex_hull_perimeter()` | Perimeter of convex hull           |
| `bounding_rect_area()`    | Area of bounding rectangle         |
| `inscribed_radius()`      | Radius of maximum inscribed circle |
| `major_axis()`            | Major axis length                  |
| `minor_axis()`            | Minor axis length                  |

### Spatial Relationships

| Function | Description |
|----|----|
| `shared_boundary()` | Total length of boundary shared with a reference polygon |

## Parallel Processing

All functions accept an `ncores` argument to enable parallel processing across polygons:

``` r
# Use 4 cores
nc <- calc_multiple_metrics(nc, ncores = 4)

# Sequential (default)
nc <- calc_multiple_metrics(nc, ncores = 1)
```

The package uses `future.apply` and automatically adjusts the number of cores if the requested value exceeds the available cores.

## Handling Geographic (Lat/Lon) Data

If your data has a geographic CRS (latitude/longitude), the package:

- Emits a warning recommending projection to a planar CRS.
- Automatically disables s2 spherical geometry for metrics that rely on `st_point_on_surface` or `st_centroid` (e.g., `boyce_clark`, `cohesion`, `exchange`, `girth`, `inscribed_radius`, `proximity`, `skew`, `sphericity`, `shared_boundary`).
- Restores the original s2 state after each function call.

For best accuracy, project your data using `st_transform()` before computing metrics:

``` r
nc_proj <- st_transform(nc, crs = 32723)  # UTM zone 23S
nc_results <- calc_multiple_metrics(nc_proj, ncores = 4)
```

## References

The metrics are based on standard landscape ecology and spatial analysis literature:

- Boyce, R. R., & Clark, W. A. V. (1964). The concept of shape in geography. *Geographical Review*, 54(4), 561-572.
- Mandelbrot, B. B. (1983). *The Fractal Geometry of Nature*. W. H. Freeman.
- McGarigal, K., Cushman, S. A., & Ene, E. (2012). *FRAGSTATS v4: Spatial Pattern Analysis Program for Categorical and Continuous Maps*. University of Massachusetts, Amherst.
- Polsby, D. D., & Popper, R. D. (1991). The third criterion: Compactness as a procedural safeguard against partisan gerrymandering. *Yale Law & Policy Review*, 9(2), 301-353.
- Reock, E. C. (1961). A note: Measuring compactness as a requirement of legislative apportionment. *Midwest Journal of Political Science*, 5(1), 70-74.
- Rosin, P. L. (2000). Measuring shape: ellipticity, rectangularity, and triangularity. *Proceedings of the 15th International Conference on Pattern Recognition*.
- Schumaker, N. H. (1996). Using landscape indices to predict habitat connectivity. *Ecology*, 77(4), 1210-1225.
- Schwartzberg, J. E. (1966). Reapportionment, gerrymanders, and the 1964 elections. *The Annals of the American Academy of Political and Social Science*, 363(1), 54-68.

## Dependencies

- `sf` – for spatial data handling.
- `future` & `future.apply` – for parallel processing.

## Contributing

Contributions are welcome! Please open an issue or submit a pull request on GitHub.

## License

This package is licensed under the MIT License. See the LICENSE file for details.

## Acknowledgements

This package was inspired by the `vectormetrics` and `redistmetrics` packages, and re-implements their metrics from scratch using base `sf` functions.

## Citation

If you use shapeMetrics in your research, please cite it as:

``` bibtex
@software{shapeMetrics2026,
  author = {Guilherme Correia},
  title = {shapeMetrics: Shape and Compactness Metrics for Spatial Polygons},
  year = {2026},
  url = {https://github.com/guilhermegcorreia9-sketch/shapeMetrics}
}
```
