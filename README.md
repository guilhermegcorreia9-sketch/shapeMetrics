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

## Function Reference

All functions below accept an `sf_obj` (an `sf` object with polygon/multipolygon geometries), an `ncores` argument for parallel processing (default `1`), and a `quiet` argument to suppress core-usage messages.

### Geometric Attributes

| Function | Description |
|----|----|
| `area()` | Area of each polygon. |
| `perimeter()` | Perimeter of each polygon. |
| `convex_hull_area()` | Area of the convex hull of each polygon. |
| `convex_hull_perimeter()` | Perimeter of the convex hull of each polygon. |
| `bounding_rect_area()` | Area of the axis-aligned minimum bounding rectangle. |
| `inscribed_radius()` | Radius of the maximum inscribed circle (approximated from boundary points). |
| `major_axis()` | Length of the major axis, based on the bounding box. |
| `minor_axis()` | Length of the minor axis, based on the bounding box. |
| `polradius()` | POLRADIUS index: maximum distance from a boundary vertex to the polygon's point-on-surface. |
| `p_gyradius()` | Average distance from all boundary vertices to the polygon's point-on-surface (radius of gyration). |
| `p_density()` | PDENSITY index: polygon area divided by `polradius`. |

### Shape Metrics (`vectormetrics` style)

| Function | Formula | Description |
|----|----|----|
| `circularity()` | A/(P²/4π) | Equivalent to the inverse of the squared shape index. |
| `compactness()` | √(4πA)/P | Form factor. Ranges 0–1, with 1 for a perfect circle. |
| `convexity()` | P_hull/P | Perimeter of the convex hull divided by the polygon perimeter. |
| `elongation()` | minor axis / major axis | Ratio of minor axis to major axis (bounding-box based). |
| `fractal_dimension()` | 2·log(P)/log(A) | Approximate fractal dimension. |
| `rectangularity()` | A/A_MABR | Polygon area divided by the area of its Minimum Area Bounding Rectangle (MABR). |
| `roughness()` | shape_index − 1 | Shape Index minus 1. |
| `shape_index()` | P/(2√(πA)) | Minimum value of 1 for a circle. |
| `convex_hull_compact()` | A/A_hull | Solidity: polygon area divided by convex hull area. |
| `sphericity()` | r_inscribed/r_circumscribed | Ratio of the maximum inscribed circle radius to the minimum circumscribing circle radius. |

### Compactness Metrics (`redistmetrics` style)

| Function | Formula | Description |
|----|----|----|
| `polsby_popper()` | 4πA/P² | Ranges 0 (least compact) to 1 (perfect circle). |
| `schwartzberg()` | 1/(P/(2π√(A/π))) | Inverse of the perimeter-to-equal-area-circle-perimeter ratio. |
| `reock()` | A/A_circumscribing circle | Polygon area divided by the area of its minimum circumscribing circle. |
| `box_reock()` | A/A_bounding rectangle | Polygon area divided by the area of its axis-aligned bounding rectangle. |
| `boyce_clark()` | Σ‖(r_i/Σr)·100 − 100/n‖ | Radial-distance-from-centroid index across `n_radii` angular sectors; lower values indicate higher compactness. |
| `skew()` | A_inscribed/A_circumscribed | Ratio of the maximum inscribed circle area to the minimum circumscribing circle area. |
| `girth()` | r_inscribed/r_equal-area | Ratio of the maximum inscribed circle radius to the equal-area circle radius. |
| `detour()` | P_equal-area circle/P_hull | Ratio of the equal-area circle perimeter to the convex hull perimeter. |

### Additional Indices

| Function | Formula | Description |
|----|----|----|
| `cohesion()` | (r²/2) / mean(d²) | Ratio of the mean squared centroid distance in an equal-area circle to that within the polygon (Monte Carlo sampled). |
| `exchange()` | A_inside equal-area circle / A | Proportion of the polygon's sampled area that falls within the equal-area circle centered on its centroid. |
| `proximity()` | (2r/3) / mean(d) | Ratio of the mean centroid distance in the equal-area circle to that within the polygon (Monte Carlo sampled). |
| `range_index()` | 2r_equal-area / 2r_circumscribed | Ratio of the equal-area circle diameter to the minimum circumscribing circle diameter. |
| `equiv_rectangular()` | P_square/P | Ratio of the perimeter of an equal-area square to the polygon's actual perimeter. |
| `perimeter_index()` | 1/shape_index | Inverse of the shape index. |
| `perimeter_area_ratio()` | P/A | Lower values indicate a more compact shape. |

### Spatial Relationships

| Function | Description |
|----|----|
| `shared_boundary(sf_obj, ref)` | Total length of each polygon's boundary shared with a single reference polygon (`ref`). |

### Batch / Group Functions

| Function | Description |
|----|----|
| `calc_multiple_metrics(sf_obj, metrics = NULL, ncores = 1, progress = FALSE, ...)` | Applies a set of metric functions to `sf_obj`, adding one column per metric. By default computes all implemented metrics; accepts a named list of functions or a character vector of metric names to compute a subset. Supports an optional progress bar. |
| `compactness_metrics(sf_obj, ncores = 1, ...)` | Convenience wrapper that computes the `redistmetrics`-style compactness subset: `polsby_popper`, `schwartzberg`, `reock`, `box_reock`, `convex_hull_compact`, `skew`, `elongation`, `boyce_clark`, `compactness`, `perimeter_index`, `detour`, `girth`, `range_index`. |
| `shape_metrics(sf_obj, ncores = 1, ...)` | Convenience wrapper that computes the `vectormetrics`-style shape subset: `polsby_popper`, `circularity`, `reock`, `box_reock`, `convex_hull_compact`, `convexity`, `sphericity`, `elongation`, `fractal_dimension`, `shape_index`, `roughness`, `compactness`, `rectangularity`, `perimeter_area_ratio`, `perimeter_index`, `cohesion`, `detour`, `equiv_rectangular`, `exchange`, `girth`, `proximity`, `range_index`. |

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
- Automatically disables s2 spherical geometry for metrics that rely on `st_point_on_surface` or `st_centroid` (e.g., `boyce_clark`, `cohesion`, `exchange`, `girth`, `inscribed_radius`, `proximity`, `skew`, `sphericity`, `shared_boundary`, `p_gyradius`, `polradius`, `p_density`).
- Restores the original s2 state after each function call.

For best accuracy, project your data using `st_transform()` before computing metrics:

``` r
nc_proj <- st_transform(nc, crs = 32723)  # UTM zone 23S
nc_results <- calc_multiple_metrics(nc_proj, ncores = 4)
```

## References

The metrics are based on standard landscape ecology and spatial analysis literature:

- Boyce, R. R., & Clark, W. A. V. (1964). The concept of shape in geography. *Geographical Review*, 54(4), 561-572.
- Körting, T. S., Fonseca, L. M. G., & Câmara, G. (2013). GeoDMA — Geographic Data Mining Analyst. *Computers & Geosciences*, 57, 133-145.
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
