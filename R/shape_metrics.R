#' @references
#' \itemize{
#'   \item{**Inspired by packages:**
#'     \itemize{
#'       \item{vectormetrics: \url{https://r-spatialecology.github.io/vectormetrics/}}
#'       \item{redistmetrics: \url{https://alarm-redist.org/redistmetrics/}}
#'     }
#'   }
#'   \item{**Scientific references for the metrics:**
#'     \itemize{
#'       \item{Boyce, R. R., & Clark, W. A. V. (1964). The concept of shape in geography. *Geographical Review*, 54(4), 561-572.}
#'       \item{Mandelbrot, B. B. (1983). *The Fractal Geometry of Nature*. W. H. Freeman.}
#'       \item{McGarigal, K., Cushman, S. A., & Ene, E. (2012). *FRAGSTATS v4: Spatial Pattern Analysis Program for Categorical and Continuous Maps*. University of Massachusetts, Amherst.}
#'       \item{Polsby, D. D., & Popper, R. D. (1991). The third criterion: Compactness as a procedural safeguard against partisan gerrymandering. *Yale Law & Policy Review*, 9(2), 301-353.}
#'       \item{Reock, E. C. (1961). A note: Measuring compactness as a requirement of legislative apportionment. *Midwest Journal of Political Science*, 5(1), 70-74.}
#'       \item{Rosin, P. L. (2000). Measuring shape: ellipticity, rectangularity, and triangularity. *Proceedings of the 15th International Conference on Pattern Recognition*.}
#'       \item{Schumaker, N. H. (1996). Using landscape indices to predict habitat connectivity. *Ecology*, 77(4), 1210-1225.}
#'       \item{Schwartzberg, J. E. (1966). Reapportionment, gerrymanders, and the 1964 elections. *The Annals of the American Academy of Political and Social Science*, 363(1), 54-68.}
#'     }
#'   }
#' }
#' @keywords internal
"_PACKAGE"

#' @import sf
#' @import future
#' @import future.apply

# Declare global variables to avoid R CMD check notes
utils::globalVariables(c("st_oriented_bbox"))

# =============================================================================
# INTERNAL GEOMETRIC HELPERS
# =============================================================================

#' @keywords internal
get_area_perimeter <- function(poly) {
  list(area = as.numeric(st_area(poly)), 
       perimeter = as.numeric(st_perimeter(poly)))
}

#' @keywords internal
get_convex_hull <- function(poly) st_convex_hull(st_union(poly))

#' @keywords internal
get_min_circumscribing_circle <- function(poly) {
  bbox <- st_bbox(poly)
  cx <- (bbox$xmin + bbox$xmax) / 2
  cy <- (bbox$ymin + bbox$ymax) / 2
  radius <- max(bbox$xmax - cx, bbox$ymax - cy)
  st_buffer(st_sfc(st_point(c(cx, cy)), crs = st_crs(poly)), dist = radius)
}

#' @keywords internal
get_min_bounding_rect <- function(poly) {
  bbox <- st_bbox(poly)
  st_sf(geometry = st_as_sfc(bbox), crs = st_crs(poly))
}

#' @keywords internal
get_mabr <- function(poly) {
  if (exists("st_oriented_bbox", where = asNamespace("sf"))) {
    st_oriented_bbox(poly)
  } else {
    warning("st_oriented_bbox not available. Using axis-aligned bounding box for rectangularity.")
    st_as_sfc(st_bbox(poly))
  }
}

#' Get a single point on surface
#'
#' Tries multiple approaches to obtain a single POINT inside the polygon:
#' 1. st_point_on_surface
#' 2. st_centroid (with of_largest_polygon = TRUE)
#' 3. st_sample (one random point)
#'
#' Returns a single POINT geometry (sfg) or NA with a warning.
#' @param poly An sf polygon geometry.
#' @return An `sfg` POINT object, or NA with a warning.
#' @keywords internal
get_point_on_surface <- function(poly) {
  # Helper to check if an object is a valid single POINT
  is_valid_point <- function(x) {
    if (is.null(x)) return(FALSE)
    if (inherits(x, "sfc")) {
      if (length(x) != 1) return(FALSE)
      x <- x[[1]]
    }
    if (!inherits(x, "sfg")) return(FALSE)
    if (!inherits(x, "POINT")) return(FALSE)
    # Ensure it has exactly 2 coordinates (XY) or 3 (XYZ) – but at least 2
    if (length(x) < 2) return(FALSE)
    # Check if coordinates are finite
    if (any(!is.finite(x))) return(FALSE)
    return(TRUE)
  }
  
  # Helper to check if point is inside polygon
  is_inside <- function(pt, poly) {
    if (is.na(pt) || !is_valid_point(pt)) return(FALSE)
    # Convert to sfc for st_contains
    if (inherits(pt, "sfg")) pt <- st_sfc(pt, crs = st_crs(poly))
    result <- tryCatch(st_contains(poly, pt, sparse = FALSE)[1, 1], error = function(e) FALSE)
    return(isTRUE(result))
  }
  
  # Try 1: st_point_on_surface (guarantees point inside)
  pts <- tryCatch(st_point_on_surface(poly), error = function(e) NULL)
  if (is_valid_point(pts) && is_inside(pts, poly)) {
    if (inherits(pts, "sfc")) pts <- pts[[1]]
    return(pts)
  }
  
  # Try 2: st_centroid (with largest polygon)
  pts <- tryCatch(st_centroid(poly, of_largest_polygon = TRUE), error = function(e) NULL)
  if (is_valid_point(pts) && is_inside(pts, poly)) {
    if (inherits(pts, "sfc")) pts <- pts[[1]]
    return(pts)
  }
  
  # Try 3: sample points inside polygon (st_sample already guarantees inside)
  pts <- tryCatch(st_sample(poly, size = 1, type = "random"), error = function(e) NULL)
  if (is_valid_point(pts) && is_inside(pts, poly)) {
    if (inherits(pts, "sfc")) pts <- pts[[1]]
    return(pts)
  }
  
  # If all fail
  warning("Could not obtain a valid point on surface. Returning NA.")
  return(NA)
}

#' @keywords internal
get_max_inscribed_circle_radius <- function(poly) {
  centroid <- get_point_on_surface(poly)
  if (!inherits(centroid, "POINT")) {
    warning("Could not compute centroid point. Returning NA.")
    return(NA)
  }
  boundary <- st_boundary(poly)
  if (length(boundary) == 0) {
    warning("Polygon has no boundary. Returning NA.")
    return(NA)
  }
  pts <- st_cast(boundary, "POINT")
  if (length(pts) < 3) {
    warning("Insufficient boundary points (<3). Returning NA.")
    return(NA)
  }
  coords <- st_coordinates(pts)
  cent <- st_coordinates(centroid)
  dists <- sqrt((coords[,1] - cent[1])^2 + (coords[,2] - cent[2])^2)
  if (length(dists) == 0 || all(is.na(dists))) {
    warning("No distances computed. Returning NA.")
    return(NA)
  }
  return(min(dists, na.rm = TRUE))
}

#' @keywords internal
get_axes <- function(poly) {
  bbox <- st_bbox(poly)
  list(major = max(bbox$xmax - bbox$xmin, bbox$ymax - bbox$ymin),
       minor = min(bbox$xmax - bbox$xmin, bbox$ymax - bbox$ymin))
}

#' @keywords internal
get_equal_area_circle_radius <- function(area) {
  sqrt(area / pi)
}

#' @keywords internal
sample_points_in_poly <- function(poly, n = 1000) {
  pts <- st_sample(poly, size = n, type = "random")
  if (length(pts) == 0) {
    warning("Sampling returned no points. Returning empty.")
    return(st_sfc())
  }
  pts
}

# =============================================================================
# SF CLASS CHECK, VALIDITY, AND S2 CONTROL
# =============================================================================

#' Check if object is an sf object and check CRS
#'
#' @param obj Object to check.
#' @param arg_name Name of the argument (for error message).
#' @return Invisibly returns TRUE if object is sf; otherwise stops with error.
#' @keywords internal
check_sf <- function(obj, arg_name = "sf_obj") {
  if (!inherits(obj, "sf")) {
    stop(sprintf("The object provided as '%s' is not an sf object. Please provide an sf object with polygon geometries.", arg_name))
  }
  crs <- st_crs(obj)
  if (is.null(crs) || is.na(crs)) {
    stop("The provided 'sf' object has an undefined or invalid CRS.")
  }
  if (st_is_longlat(crs)) {
    warning("Input data has a geographic CRS (latitude/longitude). Measurements will be computed on the ellipsoid (WGS84). For more accurate planar metrics, consider projecting to a local projected CRS using st_transform().")
  }
  invisible(TRUE)
}

#' Check if all geometries in an sf object are valid
#'
#' @param sf_obj An sf object.
#' @return Invisibly returns TRUE if all geometries are valid; otherwise stops with error.
#' @keywords internal
check_valid_geometry <- function(sf_obj) {
  invalid <- !st_is_valid(sf_obj)
  if (any(invalid)) {
    invalid_indices <- which(invalid)
    stop(sprintf(
      "Invalid geometries detected in polygons: %s. Please fix invalid geometries before computing metrics.",
      paste(invalid_indices, collapse = ", ")
    ))
  }
  invisible(TRUE)
}

#' Check if all geometries in an sf object are polygons
#'
#' @param sf_obj An sf object.
#' @param type Valid geometry type.
#' @return Invisibly returns TRUE if all geometries are polygons; otherwise stops with error.
#' @keywords internal
check_polygon_geometry <- function(sf_obj, type = "POLYGON") {
  is_polygon <- sf::st_is(sf_obj, type)
  
  if (any(!is_polygon)) {
    invalid_indices <- which(!is_polygon)
    type_formatted <- paste(type, collapse = " or ")
    stop(sprintf(
      "Non-polygon geometries detected at indices: %s. Please provide only %s geometries.",
      paste(invalid_indices, collapse = ", "),
      type_formatted
    ))
  }
  
  invisible(TRUE)
}

#' Disable s2 for geographic CRS and return previous state
#'
#' If the sf object has a geographic CRS and s2 is currently enabled,
#' this function disables s2 and informs the user. It returns the previous
#' state of s2 so it can be restored later.
#'
#' @param sf_obj An sf object.
#' @param func_name Name of the calling function (for the message).
#' @return A logical value: the previous state of s2 (TRUE if s2 was enabled), or NULL if no change.
#' @keywords internal
disable_s2_if_geographic <- function(sf_obj, func_name = "function") {
  check_sf(sf_obj)
  crs <- st_crs(sf_obj)
  if (is.na(crs) || is.null(crs)) {
    return(invisible(NULL))
  }
  if (st_is_longlat(crs)) {
    previous_state <- sf_use_s2()
    if (previous_state) {
      sf_use_s2(FALSE)
      message(sprintf(
        "Geographic CRS detected in %s. Disabling spherical geometry (s2) to avoid issues with point-on-surface calculations. For best accuracy, consider projecting your data to a local projected CRS using st_transform().",
        func_name
      ))
    }
    return(invisible(previous_state))
  }
  return(invisible(NULL))
}

# =============================================================================
# SCALAR FUNCTIONS (internal) – in alphabetical order
# =============================================================================

#' @keywords internal
calc_area <- function(poly) {
  val <- as.numeric(st_area(poly))
  if (is.na(val) || val == 0) {
    warning("Area is zero or NA. Returning NA.")
    return(NA)
  }
  val
}

#' @keywords internal
calc_bounding_rect_area <- function(poly) {
  rect <- get_min_bounding_rect(poly)
  if (is.na(rect)) {
    warning("Could not compute bounding rectangle. Returning NA.")
    return(NA)
  }
  val <- as.numeric(st_area(rect))
  if (is.na(val) || val == 0) {
    warning("Bounding rectangle area is zero or NA. Returning NA.")
    return(NA)
  }
  val
}

#' @keywords internal
calc_box_reock <- function(poly) {
  area <- as.numeric(st_area(poly))
  if (area == 0 || is.na(area)) {
    warning("Zero or NA area. Returning NA.")
    return(NA)
  }
  rect <- get_min_bounding_rect(poly)
  if (is.na(rect)) {
    warning("Could not compute bounding rectangle. Returning NA.")
    return(NA)
  }
  area_rect <- as.numeric(st_area(rect))
  if (area_rect == 0 || is.na(area_rect)) {
    warning("Bounding rectangle area is zero or NA. Returning NA.")
    return(NA)
  }
  area / area_rect
}

#' @keywords internal
calc_boyce_clark <- function(poly, n_radii = 16) {
  centroid <- get_point_on_surface(poly)
  if (!inherits(centroid, "POINT")) {
    warning("Could not obtain a valid point on surface. Returning NA.")
    return(NA)
  }
  boundary <- st_boundary(poly)
  if (length(boundary) == 0) {
    warning("Polygon has no boundary. Returning NA.")
    return(NA)
  }
  pts <- st_cast(boundary, "POINT")
  if (length(pts) < 3) {
    warning("Insufficient boundary points (<3). Returning NA.")
    return(NA)
  }
  coords <- st_coordinates(pts)
  cent <- st_coordinates(centroid)
  dx <- coords[,1] - cent[1]; dy <- coords[,2] - cent[2]
  dists <- sqrt(dx^2 + dy^2)
  angles <- atan2(dy, dx)
  breaks <- seq(-pi, pi, length.out = n_radii + 1)
  groups <- cut(angles, breaks, include.lowest = TRUE)
  radii <- tapply(dists, groups, mean, na.rm = TRUE)
  radii[is.na(radii)] <- 0
  total <- sum(radii)
  if (total == 0 || is.na(total)) {
    warning("Total radial sum is zero or NA. Returning NA.")
    return(NA)
  }
  expected <- 100 / n_radii
  sum(abs((radii / total) * 100 - expected))
}

#' @keywords internal
calc_circularity <- function(poly) {
  ap <- get_area_perimeter(poly)
  if (ap$area == 0 || ap$perimeter == 0 || is.na(ap$area) || is.na(ap$perimeter)) {
    warning("Zero or NA area or perimeter. Returning NA.")
    return(NA)
  }
  ap$area / ((ap$perimeter^2) / (4 * pi))
}

#' @keywords internal
calc_cohesion <- function(poly, n_samples = 1000) {
  area <- as.numeric(st_area(poly))
  if (area == 0 || is.na(area)) {
    warning("Zero or NA area. Returning NA.")
    return(NA)
  }
  r <- get_equal_area_circle_radius(area)
  pts <- sample_points_in_poly(poly, n = n_samples)
  if (length(pts) == 0) {
    warning("Sampling returned no points. Returning NA.")
    return(NA)
  }
  centroid <- get_point_on_surface(poly)
  if (!inherits(centroid, "POINT")) {
    warning("Could not obtain a valid point on surface. Returning NA.")
    return(NA)
  }
  coords <- st_coordinates(pts)
  cent <- st_coordinates(centroid)
  if (any(is.na(cent))) {
    warning("Could not compute centroid coordinates. Returning NA.")
    return(NA)
  }
  d_sq_poly <- (coords[,1] - cent[1])^2 + (coords[,2] - cent[2])^2
  mean_d_sq_poly <- mean(d_sq_poly, na.rm = TRUE)
  if (is.na(mean_d_sq_poly) || mean_d_sq_poly == 0) {
    warning("Mean squared distance is zero or NA. Returning NA.")
    return(NA)
  }
  mean_d_sq_circle <- r^2 / 2
  mean_d_sq_circle / mean_d_sq_poly
}

#' @keywords internal
calc_compactness <- function(poly) {
  ap <- get_area_perimeter(poly)
  if (ap$area == 0 || ap$perimeter == 0 || is.na(ap$area) || is.na(ap$perimeter)) {
    warning("Zero or NA area or perimeter. Returning NA.")
    return(NA)
  }
  sqrt(4 * pi * ap$area) / ap$perimeter
}

#' @keywords internal
calc_convex_hull_area <- function(poly) {
  hull <- get_convex_hull(poly)
  if (is.na(hull) || length(hull) == 0) {
    warning("Could not compute convex hull. Returning NA.")
    return(NA)
  }
  val <- as.numeric(st_area(hull))
  if (is.na(val) || val == 0) {
    warning("Convex hull area is zero or NA. Returning NA.")
    return(NA)
  }
  val
}

#' @keywords internal
calc_convex_hull_compact <- function(poly) {
  area <- as.numeric(st_area(poly))
  if (area == 0 || is.na(area)) {
    warning("Zero or NA area. Returning NA.")
    return(NA)
  }
  hull <- get_convex_hull(poly)
  if (is.na(hull) || length(hull) == 0) {
    warning("Could not compute convex hull. Returning NA.")
    return(NA)
  }
  area_hull <- as.numeric(st_area(hull))
  if (area_hull == 0 || is.na(area_hull)) {
    warning("Convex hull area is zero or NA. Returning NA.")
    return(NA)
  }
  area / area_hull
}

#' @keywords internal
calc_convex_hull_perimeter <- function(poly) {
  hull <- get_convex_hull(poly)
  if (is.na(hull) || length(hull) == 0) {
    warning("Could not compute convex hull. Returning NA.")
    return(NA)
  }
  val <- as.numeric(st_perimeter(hull))
  if (is.na(val) || val == 0) {
    warning("Convex hull perimeter is zero or NA. Returning NA.")
    return(NA)
  }
  val
}

#' @keywords internal
calc_convexity <- function(poly) {
  ap <- get_area_perimeter(poly)
  if (ap$perimeter == 0 || is.na(ap$perimeter)) {
    warning("Zero or NA perimeter. Returning NA.")
    return(NA)
  }
  hull <- get_convex_hull(poly)
  if (is.na(hull) || length(hull) == 0) {
    warning("Could not compute convex hull. Returning NA.")
    return(NA)
  }
  perim_hull <- as.numeric(st_perimeter(hull))
  if (perim_hull == 0 || is.na(perim_hull)) {
    warning("Convex hull perimeter is zero or NA. Returning NA.")
    return(NA)
  }
  perim_hull / ap$perimeter
}

#' @keywords internal
calc_detour <- function(poly) {
  area <- as.numeric(st_area(poly))
  if (area == 0 || is.na(area)) {
    warning("Zero or NA area. Returning NA.")
    return(NA)
  }
  hull <- get_convex_hull(poly)
  if (is.na(hull) || length(hull) == 0) {
    warning("Could not compute convex hull. Returning NA.")
    return(NA)
  }
  perim_hull <- as.numeric(st_perimeter(hull))
  if (perim_hull == 0 || is.na(perim_hull)) {
    warning("Convex hull perimeter is zero or NA. Returning NA.")
    return(NA)
  }
  perim_circle <- 2 * pi * sqrt(area / pi)
  perim_circle / perim_hull
}

#' @keywords internal
calc_elongation <- function(poly) {
  axes <- get_axes(poly)
  if (axes$major == 0 || is.na(axes$major)) {
    warning("Major axis is zero or NA. Returning NA.")
    return(NA)
  }
  axes$minor / axes$major
}

#' @keywords internal
calc_equiv_rectangular <- function(poly) {
  ap <- get_area_perimeter(poly)
  if (ap$area == 0 || ap$perimeter == 0 || is.na(ap$area) || is.na(ap$perimeter)) {
    warning("Zero or NA area or perimeter. Returning NA.")
    return(NA)
  }
  side <- sqrt(ap$area)
  perim_square <- 4 * side
  perim_square / ap$perimeter
}

#' @keywords internal
calc_exchange <- function(poly, n_samples = 1000) {
  area <- as.numeric(st_area(poly))
  if (area == 0 || is.na(area)) {
    warning("Zero or NA area. Returning NA.")
    return(NA)
  }
  r <- sqrt(area / pi)
  centroid <- get_point_on_surface(poly)
  if (!inherits(centroid, "POINT")) {
    warning("Could not obtain a valid point on surface. Returning NA.")
    return(NA)
  }
  pts <- st_sample(poly, size = n_samples, type = "random")
  if (length(pts) == 0) {
    warning("Sampling returned no points. Returning NA.")
    return(NA)
  }
  coords <- st_coordinates(pts)
  cent <- st_coordinates(centroid)
  if (any(is.na(cent))) {
    warning("Could not compute centroid coordinates. Returning NA.")
    return(NA)
  }
  dists <- sqrt((coords[,1] - cent[1])^2 + (coords[,2] - cent[2])^2)
  inside <- sum(dists <= r, na.rm = TRUE)
  if (length(pts) == 0) {
    warning("No points sampled. Returning NA.")
    return(NA)
  }
  inside / length(pts)
}

#' @keywords internal
calc_fractal_dimension <- function(poly) {
  ap <- get_area_perimeter(poly)
  if (ap$area == 0 || ap$perimeter == 0 || is.na(ap$area) || is.na(ap$perimeter)) {
    warning("Zero or NA area or perimeter. Returning NA.")
    return(NA)
  }
  2 * log(ap$perimeter) / log(ap$area)
}

#' @keywords internal
calc_girth <- function(poly) {
  area <- as.numeric(st_area(poly))
  if (area == 0 || is.na(area)) {
    warning("Zero or NA area. Returning NA.")
    return(NA)
  }
  r_inscribed <- get_max_inscribed_circle_radius(poly)
  if (is.na(r_inscribed) || r_inscribed == 0) {
    warning("Maximum inscribed radius is zero or NA. Returning NA.")
    return(NA)
  }
  r_equal <- sqrt(area / pi)
  if (r_equal == 0 || is.na(r_equal)) {
    warning("Equal-area radius is zero or NA. Returning NA.")
    return(NA)
  }
  r_inscribed / r_equal
}

#' @keywords internal
calc_inscribed_radius <- function(poly) {
  get_max_inscribed_circle_radius(poly)
}

#' @keywords internal
calc_major_axis <- function(poly) {
  get_axes(poly)$major
}

#' @keywords internal
calc_minor_axis <- function(poly) {
  get_axes(poly)$minor
}

#' @keywords internal
calc_perimeter <- function(poly) {
  val <- as.numeric(st_perimeter(poly))
  if (is.na(val) || val == 0) {
    warning("Perimeter is zero or NA. Returning NA.")
    return(NA)
  }
  val
}

#' @keywords internal
calc_perimeter_area_ratio <- function(poly) {
  ap <- get_area_perimeter(poly)
  if (ap$area == 0 || is.na(ap$area)) {
    warning("Zero or NA area. Returning NA.")
    return(NA)
  }
  ap$perimeter / ap$area
}

#' @keywords internal
calc_perimeter_index <- function(poly) {
  si <- calc_shape_index(poly)
  if (is.na(si) || si == 0) {
    warning("Shape index is zero or NA. Returning NA.")
    return(NA)
  }
  1 / si
}

#' @keywords internal
calc_polsby_popper <- function(poly) {
  ap <- get_area_perimeter(poly)
  if (ap$perimeter == 0 || is.na(ap$perimeter)) {
    warning("Zero or NA perimeter. Returning NA.")
    return(NA)
  }
  (4 * pi * ap$area) / (ap$perimeter^2)
}

#' @keywords internal
calc_proximity <- function(poly, n_samples = 1000) {
  area <- as.numeric(st_area(poly))
  if (area == 0 || is.na(area)) {
    warning("Zero or NA area. Returning NA.")
    return(NA)
  }
  r <- sqrt(area / pi)
  centroid <- get_point_on_surface(poly)
  if (!inherits(centroid, "POINT")) {
    warning("Could not obtain a valid point on surface. Returning NA.")
    return(NA)
  }
  pts <- st_sample(poly, size = n_samples, type = "random")
  if (length(pts) == 0) {
    warning("Sampling returned no points. Returning NA.")
    return(NA)
  }
  coords <- st_coordinates(pts)
  cent <- st_coordinates(centroid)
  if (any(is.na(cent))) {
    warning("Could not compute centroid coordinates. Returning NA.")
    return(NA)
  }
  dists <- sqrt((coords[,1] - cent[1])^2 + (coords[,2] - cent[2])^2)
  mean_dist_poly <- mean(dists, na.rm = TRUE)
  if (is.na(mean_dist_poly) || mean_dist_poly == 0) {
    warning("Mean distance is zero or NA. Returning NA.")
    return(NA)
  }
  mean_dist_circle <- 2 * r / 3
  mean_dist_circle / mean_dist_poly
}

#' @keywords internal
calc_range_index <- function(poly) {
  area <- as.numeric(st_area(poly))
  if (area == 0 || is.na(area)) {
    warning("Zero or NA area. Returning NA.")
    return(NA)
  }
  circle <- get_min_circumscribing_circle(poly)
  if (is.na(circle) || length(circle) == 0) {
    warning("Could not compute circumscribing circle. Returning NA.")
    return(NA)
  }
  bbox <- st_bbox(circle)
  r_out <- max(bbox$xmax - bbox$xmin, bbox$ymax - bbox$ymin) / 2
  if (r_out == 0 || is.na(r_out)) {
    warning("Circumscribing radius is zero or NA. Returning NA.")
    return(NA)
  }
  r_equal <- sqrt(area / pi)
  (2 * r_equal) / (2 * r_out)
}

#' @keywords internal
calc_rectangularity <- function(poly) {
  area <- as.numeric(st_area(poly))
  if (area == 0 || is.na(area)) {
    warning("Zero or NA area. Returning NA.")
    return(NA)
  }
  mabr <- get_mabr(poly)
  if (is.na(mabr) || length(mabr) == 0) {
    warning("Could not compute MABR. Returning NA.")
    return(NA)
  }
  area_mabr <- as.numeric(st_area(mabr))
  if (area_mabr == 0 || is.na(area_mabr)) {
    warning("MABR area is zero or NA. Returning NA.")
    return(NA)
  }
  area / area_mabr
}

#' @keywords internal
calc_reock <- function(poly) {
  area <- as.numeric(st_area(poly))
  if (area == 0 || is.na(area)) {
    warning("Zero or NA area. Returning NA.")
    return(NA)
  }
  circle <- get_min_circumscribing_circle(poly)
  if (is.na(circle) || length(circle) == 0) {
    warning("Could not compute circumscribing circle. Returning NA.")
    return(NA)
  }
  area_circle <- as.numeric(st_area(circle))
  if (area_circle == 0 || is.na(area_circle)) {
    warning("Circumscribing circle area is zero or NA. Returning NA.")
    return(NA)
  }
  area / area_circle
}

#' @keywords internal
calc_roughness <- function(poly) {
  si <- calc_shape_index(poly)
  if (is.na(si)) {
    warning("Shape index is NA. Returning NA.")
    return(NA)
  }
  si - 1
}

#' @keywords internal
calc_schwartzberg <- function(poly) {
  ap <- get_area_perimeter(poly)
  if (ap$area == 0 || is.na(ap$area)) {
    warning("Zero or NA area. Returning NA.")
    return(NA)
  }
  denom <- ap$perimeter / (2 * pi * sqrt(ap$area / pi))
  if (denom == 0 || is.na(denom)) {
    warning("Denominator is zero or NA. Returning NA.")
    return(NA)
  }
  1 / denom
}

#' @keywords internal
calc_shape_index <- function(poly) {
  ap <- get_area_perimeter(poly)
  if (ap$area == 0 || is.na(ap$area)) {
    warning("Zero or NA area. Returning NA.")
    return(NA)
  }
  ap$perimeter / (2 * sqrt(pi * ap$area))
}

#' @keywords internal
calc_skew <- function(poly) {
  r_in <- get_max_inscribed_circle_radius(poly)
  if (is.na(r_in) || r_in == 0) {
    warning("Maximum inscribed radius is zero or NA. Returning NA.")
    return(NA)
  }
  circle <- get_min_circumscribing_circle(poly)
  if (is.na(circle) || length(circle) == 0) {
    warning("Could not compute circumscribing circle. Returning NA.")
    return(NA)
  }
  bbox <- st_bbox(circle)
  r_out <- max(bbox$xmax - bbox$xmin, bbox$ymax - bbox$ymin) / 2
  if (r_out == 0 || is.na(r_out)) {
    warning("Circumscribing radius is zero or NA. Returning NA.")
    return(NA)
  }
  (r_in^2) / (r_out^2)
}

#' @keywords internal
calc_sphericity <- function(poly) {
  r_in <- get_max_inscribed_circle_radius(poly)
  if (is.na(r_in) || r_in == 0) {
    warning("Maximum inscribed radius is zero or NA. Returning NA.")
    return(NA)
  }
  circle <- get_min_circumscribing_circle(poly)
  if (is.na(circle) || length(circle) == 0) {
    warning("Could not compute circumscribing circle. Returning NA.")
    return(NA)
  }
  bbox <- st_bbox(circle)
  r_out <- max(bbox$xmax - bbox$xmin, bbox$ymax - bbox$ymin) / 2
  if (r_out == 0 || is.na(r_out)) {
    warning("Circumscribing radius is zero or NA. Returning NA.")
    return(NA)
  }
  r_in / r_out
}

#' @keywords internal
calc_p_gyradius <- function(poly) {
  centroid <- get_point_on_surface(poly)
  # Verifica se é um POINT válido (NÃO use is.na() em geometrias!)
  if (!inherits(centroid, "POINT")) {
    warning("Could not obtain a valid point on surface. Returning NA.")
    return(NA)
  }
  boundary <- st_boundary(poly)
  if (length(boundary) == 0) {
    warning("Polygon has no boundary. Returning NA.")
    return(NA)
  }
  pts <- st_cast(boundary, "POINT")
  if (length(pts) < 1) {
    warning("No boundary points. Returning NA.")
    return(NA)
  }
  coords <- st_coordinates(pts)
  cent <- st_coordinates(centroid)
  dists <- sqrt((coords[,1] - cent[1])^2 + (coords[,2] - cent[2])^2)
  mean(dists, na.rm = TRUE)
}

#' @keywords internal
calc_polradius <- function(poly) {
  centroid <- get_point_on_surface(poly)
  # Verifica se é um POINT válido (NÃO use is.na() em geometrias!)
  if (!inherits(centroid, "POINT")) {
    warning("Could not obtain a valid point on surface. Returning NA.")
    return(NA)
  }
  boundary <- st_boundary(poly)
  if (length(boundary) == 0) {
    warning("Polygon has no boundary. Returning NA.")
    return(NA)
  }
  pts <- st_cast(boundary, "POINT")
  if (length(pts) < 1) {
    warning("No boundary points. Returning NA.")
    return(NA)
  }
  coords <- st_coordinates(pts)
  cent <- st_coordinates(centroid)
  dists <- sqrt((coords[,1] - cent[1])^2 + (coords[,2] - cent[2])^2)
  max(dists, na.rm = TRUE)
}

#' @keywords internal
calc_p_density <- function(poly) {
  area <- as.numeric(st_area(poly))
  if (area == 0 || is.na(area)) {
    warning("Zero or NA area. Returning NA.")
    return(NA)
  }
  radius <- calc_polradius(poly)
  if (is.na(radius) || radius == 0) {
    warning("Polygon radius is zero or NA. Returning NA.")
    return(NA)
  }
  area / radius
}

# =============================================================================
# PARALLEL APPLY HELPER (com future.seed = TRUE)
# =============================================================================

#' Apply a scalar metric in parallel using future_lapply
#'
#' @param sf_obj An sf object with multiple polygons.
#' @param calc_func A scalar function (calc_*) to apply to each polygon.
#' @param col_name Name of the new column to add.
#' @param ncores Number of cores for parallel processing. If \code{>1},
#'   sets \code{future::plan("multisession", workers = ncores)}. 
#'   If \code{NULL} or \code{<1}, it is set to 1 (no parallelization).
#' @param quiet Logical: if TRUE, suppress messages about core usage.
#' @param ... Additional arguments passed to \code{calc_func}.
#' @return The sf object with the new column added.
#' @keywords internal
add_metric <- function(sf_obj, calc_func, col_name, ncores = 1, quiet = FALSE, ...) {
  check_sf(sf_obj)
  check_valid_geometry(sf_obj)
  check_polygon_geometry(sf_obj)
  
  if (is.null(ncores) || ncores < 1) ncores <- 1
  ncores <- as.integer(ncores)
  max_cores <- future::availableCores()
  if (ncores > max_cores) {
    warning(sprintf("Requested %d cores, but only %d are available. Using %d cores.", 
                    ncores, max_cores, max_cores - 1))
    ncores <- max(1, max_cores - 1)
  }
  
  if (ncores > 1 && !quiet) {
    message(sprintf("Using parallel processing with %d cores.", ncores))
  }
  
  if (ncores > 1) {
    future::plan("multisession", workers = ncores)
  }
  
  n_rows <- nrow(sf_obj)
  
  # Chunk processing thereshold
  CHUNK_THRESHOLD <- 100000
  
  values <- if (n_rows > CHUNK_THRESHOLD) {
    # Large count polygons
    chunk_size <- 5000
    n_chunks <- ceiling(n_rows / chunk_size)
  
    if (!quiet) {
      message(sprintf("Large dataset (%d polygons). Processing in %d chunks.", 
                      n_rows, n_chunks))
    }
    
    # Chunk processing function
    process_chunk <- function(chunk_idx) {
      start <- (chunk_idx - 1) * chunk_size + 1
      end <- min(chunk_idx * chunk_size, n_rows)
      chunk <- sf_obj[start:end, ]
      
      gc()
      
      # Sequential polygon processing inside the chunk
      chunk_results <- lapply(seq_len(nrow(chunk)), function(j) {
        tryCatch({
          calc_func(chunk[j, ], ...)
        }, error = function(e) {
          warning(sprintf("Error calculating metric for polygon %d: %s", 
                          start + j - 1, e$message))
          return(NA)
        })
      })
      
      # Liberar chunk da memória explicitamente
      rm(chunk)
      gc()
      
      return(chunk_results)
    }
    
    # Parallel chunk processing
    chunk_results_list <- tryCatch({
      future.apply::future_lapply(seq_len(n_chunks), process_chunk, 
                                  future.seed = TRUE, future.packages = "shapeMetrics")
    }, interrupt = function(e) {
      message("Processing interrupted by user.")
      stop(e)
    })
    
    unlist(chunk_results_list, recursive = FALSE)
  } else {
    # Low polygon count method
    tryCatch({
      future.apply::future_lapply(seq_len(n_rows), function(i) {
        tryCatch({
          calc_func(sf_obj[i, ], ...)
        }, error = function(e) {
          warning(sprintf("Error calculating metric for polygon %d: %s", i, e$message))
          return(NA)
        })
      }, future.seed = TRUE, future.packages = "shapeMetrics")
    }, interrupt = function(e) {
      message("Processing interrupted by user.")
      stop(e)
    })
  }
  
  gc()
  
  sf_obj[[col_name]] <- unlist(values)
  return(sf_obj)
}

# =============================================================================
# EXPORTED FUNCTIONS – GEOMETRIC ATTRIBUTES (alphabetical)
# =============================================================================

#' Area of polygons
#'
#' Computes the area of each polygon in the sf object.
#'
#' @param sf_obj An sf object containing polygon geometries.
#' @param ncores Number of cores for parallel computation. Default is 1 (sequential).
#' @param quiet Logical: suppress messages about core usage.
#' @return The sf object with an additional column `area`.
#' @references
#' Standard geometric measure.
#' @export
area <- function(sf_obj, ncores = 1, quiet = FALSE) {
  add_metric(sf_obj, calc_area, "area", ncores, quiet = quiet)
}

#' Area of minimum bounding rectangle
#'
#' Computes the area of the axis‑aligned minimum bounding rectangle of each polygon.
#'
#' @inheritParams area
#' @return The sf object with an additional column `bounding_rect_area`.
#' @references
#' Standard geometric measure; often used in landscape ecology.
#' @export
bounding_rect_area <- function(sf_obj, ncores = 1, quiet = FALSE) {
  add_metric(sf_obj, calc_bounding_rect_area, "bounding_rect_area", ncores, quiet = quiet)
}

#' Box‑Reock Compactness
#'
#' Area of polygon divided by area of its axis‑aligned minimum bounding rectangle.
#'
#' @inheritParams area
#' @return sf object with column `box_reock`.
#' @references
#' Reock, E. C. (1961). A note: Measuring compactness as a requirement of
#' legislative apportionment. *Midwest Journal of Political Science*, 5(1), 70-74.
#' (Axis-aligned version is a common variant.)
#' @export
box_reock <- function(sf_obj, ncores = 1, quiet = FALSE) {
  add_metric(sf_obj, calc_box_reock, "box_reock", ncores, quiet = quiet)
}

#' Boyce‑Clark Compactness
#'
#' Computes the Boyce‑Clark index based on radial distances from the centroid.
#' Lower values indicate higher compactness.
#'
#' @inheritParams area
#' @param n_radii Number of radii (rays) to use (default = 16).
#' @return The sf object with an additional column `boyce_clark`.
#' @references
#' Boyce, R. R., & Clark, W. A. V. (1964). The concept of shape in geography.
#' *Geographical Review*, 54(4), 561-572.
#' @export
boyce_clark <- function(sf_obj, ncores = 1, n_radii = 16, quiet = FALSE) {
  previous_s2 <- disable_s2_if_geographic(sf_obj, "boyce_clark")
  result <- add_metric(sf_obj, calc_boyce_clark, "boyce_clark", ncores, quiet = quiet, n_radii = n_radii)
  if (!is.null(previous_s2)) {
    sf_use_s2(previous_s2)
  }
  return(result)
}

#' Circularity
#'
#' Computes circularity: \eqn{A / (P^2 / 4\pi)}.
#' Equivalent to the inverse of the shape index squared.
#'
#' @inheritParams area
#' @return sf object with column `circularity`.
#' @references
#' McGarigal, K., Cushman, S. A., & Ene, E. (2012). FRAGSTATS v4: Spatial Pattern Analysis Program.
#' @export
circularity <- function(sf_obj, ncores = 1, quiet = FALSE) {
  add_metric(sf_obj, calc_circularity, "circularity", ncores, quiet = quiet)
}

#' Cohesion Index
#'
#' Computes the cohesion index as the ratio of the mean squared distance from
#' the centroid in an equal-area circle to that in the polygon. Higher values
#' indicate a more cohesive shape.
#'
#' @inheritParams area
#' @param n_samples Number of random points to sample inside the polygon (default: 1000).
#' @return The sf object with an additional column `cohesion`.
#' @references
#' Schumaker, N. H. (1996). Using landscape indices to predict habitat connectivity.
#' *Ecology*, 77(4), 1210-1225.
#' @export
cohesion <- function(sf_obj, ncores = 1, n_samples = 1000, quiet = FALSE) {
  previous_s2 <- disable_s2_if_geographic(sf_obj, "cohesion")
  result <- add_metric(sf_obj, calc_cohesion, "cohesion", ncores, quiet = quiet, n_samples = n_samples)
  if (!is.null(previous_s2)) {
    sf_use_s2(previous_s2)
  }
  return(result)
}

#' Compactness (Form Factor)
#'
#' Computes the compactness index (also known as form factor):
#' \deqn{\frac{\sqrt{4 \cdot \pi \cdot A}}{P}}
#' where \eqn{A} is area and \eqn{P} is perimeter. Values range from 0 to 1,
#' with 1 representing a perfect circle.
#'
#' @inheritParams area
#' @return The sf object with an additional column `compactness`.
#' @references
#' McGarigal, K., Cushman, S. A., & Ene, E. (2012). FRAGSTATS v4.
#' @export
compactness <- function(sf_obj, ncores = 1, quiet = FALSE) {
  add_metric(sf_obj, calc_compactness, "compactness", ncores, quiet = quiet)
}

#' Area of convex hull
#'
#' Computes the area of the convex hull of each polygon.
#'
#' @inheritParams area
#' @return The sf object with an additional column `convex_hull_area`.
#' @references
#' Standard geometric measure.
#' @export
convex_hull_area <- function(sf_obj, ncores = 1, quiet = FALSE) {
  add_metric(sf_obj, calc_convex_hull_area, "convex_hull_area", ncores, quiet = quiet)
}

#' Convex Hull Compactness
#'
#' Area of polygon divided by area of its convex hull.
#'
#' @inheritParams area
#' @return sf object with column `convex_hull_compact`.
#' @references
#' McGarigal, K., Cushman, S. A., & Ene, E. (2012). FRAGSTATS v4.
#' @export
convex_hull_compact <- function(sf_obj, ncores = 1, quiet = FALSE) {
  add_metric(sf_obj, calc_convex_hull_compact, "convex_hull_compact", ncores, quiet = quiet)
}

#' Perimeter of convex hull
#'
#' Computes the perimeter of the convex hull of each polygon.
#'
#' @inheritParams area
#' @return The sf object with an additional column `convex_hull_perimeter`.
#' @references
#' Standard geometric measure.
#' @export
convex_hull_perimeter <- function(sf_obj, ncores = 1, quiet = FALSE) {
  add_metric(sf_obj, calc_convex_hull_perimeter, "convex_hull_perimeter", ncores, quiet = quiet)
}

#' Convexity
#'
#' Perimeter of convex hull divided by perimeter of polygon.
#'
#' @inheritParams area
#' @return sf object with column `convexity`.
#' @references
#' McGarigal, K., Cushman, S. A., & Ene, E. (2012). FRAGSTATS v4.
#' @export
convexity <- function(sf_obj, ncores = 1, quiet = FALSE) {
  add_metric(sf_obj, calc_convexity, "convexity", ncores, quiet = quiet)
}

#' Detour Index
#'
#' Ratio of the perimeter of an equal-area circle to the perimeter of the
#' convex hull. Values close to 1 indicate a shape that is both compact and
#' convex.
#'
#' @inheritParams area
#' @return The sf object with an additional column `detour`.
#' @references
#' McGarigal, K., Cushman, S. A., & Ene, E. (2012). FRAGSTATS v4.
#' @export
detour <- function(sf_obj, ncores = 1, quiet = FALSE) {
  add_metric(sf_obj, calc_detour, "detour", ncores, quiet = quiet)
}

#' Elongation
#'
#' Ratio of minor axis to major axis (based on axis‑aligned bounding box).
#'
#' @inheritParams area
#' @return sf object with column `elongation`.
#' @references
#' McGarigal, K., Cushman, S. A., & Ene, E. (2012). FRAGSTATS v4.
#' @export
elongation <- function(sf_obj, ncores = 1, quiet = FALSE) {
  add_metric(sf_obj, calc_elongation, "elongation", ncores, quiet = quiet)
}

#' Equivalent Rectangular Index
#'
#' Ratio of the perimeter of a square with the same area as the polygon to the
#' actual perimeter. Values range from 0 to 1, with 1 representing a perfect
#' square.
#'
#' @inheritParams area
#' @return The sf object with an additional column `equiv_rectangular`.
#' @references
#' McGarigal, K., Cushman, S. A., & Ene, E. (2012). FRAGSTATS v4.
#' @export
equiv_rectangular <- function(sf_obj, ncores = 1, quiet = FALSE) {
  add_metric(sf_obj, calc_equiv_rectangular, "equiv_rectangular", ncores, quiet = quiet)
}

#' Exchange Index
#'
#' Computes the proportion of the polygon's area that lies within the equal‑area
#' circle centered on the polygon's centroid. Values close to 1 indicate a
#' shape that is well‑centered.
#'
#' @inheritParams area
#' @param n_samples Number of random points to sample inside the polygon (default: 1000).
#' @return The sf object with an additional column `exchange`.
#' @references
#' McGarigal, K., Cushman, S. A., & Ene, E. (2012). FRAGSTATS v4.
#' @export
exchange <- function(sf_obj, ncores = 1, n_samples = 1000, quiet = FALSE) {
  previous_s2 <- disable_s2_if_geographic(sf_obj, "exchange")
  result <- add_metric(sf_obj, calc_exchange, "exchange", ncores, quiet = quiet, n_samples = n_samples)
  if (!is.null(previous_s2)) {
    sf_use_s2(previous_s2)
  }
  return(result)
}

#' Fractal Dimension
#'
#' Approximate fractal dimension: \eqn{2 * \log(P) / \log(A)}.
#'
#' @inheritParams area
#' @return sf object with column `fractal_dimension`.
#' @references
#' Mandelbrot, B. B. (1983). *The Fractal Geometry of Nature*. W. H. Freeman.
#' Also used in FRAGSTATS.
#' @export
fractal_dimension <- function(sf_obj, ncores = 1, quiet = FALSE) {
  add_metric(sf_obj, calc_fractal_dimension, "fractal_dimension", ncores, quiet = quiet)
}

#' Girth Index
#'
#' Ratio of the radius of the maximum inscribed circle to the radius of the
#' equal‑area circle. Values close to 1 indicate a shape that is well‑filled.
#'
#' @inheritParams area
#' @return The sf object with an additional column `girth`.
#' @references
#' McGarigal, K., Cushman, S. A., & Ene, E. (2012). FRAGSTATS v4.
#' @export
girth <- function(sf_obj, ncores = 1, quiet = FALSE) {
  previous_s2 <- disable_s2_if_geographic(sf_obj, "girth")
  result <- add_metric(sf_obj, calc_girth, "girth", ncores, quiet = quiet)
  if (!is.null(previous_s2)) {
    sf_use_s2(previous_s2)
  }
  return(result)
}

#' Radius of maximum inscribed circle
#'
#' Approximates the radius of the maximum inscribed circle for each polygon.
#'
#' @inheritParams area
#' @return The sf object with an additional column `inscribed_radius`.
#' @references
#' Standard geometric measure.
#' @export
inscribed_radius <- function(sf_obj, ncores = 1, quiet = FALSE) {
  previous_s2 <- disable_s2_if_geographic(sf_obj, "inscribed_radius")
  result <- add_metric(sf_obj, calc_inscribed_radius, "inscribed_radius", ncores, quiet = quiet)
  if (!is.null(previous_s2)) {
    sf_use_s2(previous_s2)
  }
  return(result)
}

#' Major axis length
#'
#' Computes the length of the major axis (based on the bounding box) for each polygon.
#'
#' @inheritParams area
#' @return The sf object with an additional column `major_axis`.
#' @references
#' Standard geometric measure.
#' @export
major_axis <- function(sf_obj, ncores = 1, quiet = FALSE) {
  add_metric(sf_obj, calc_major_axis, "major_axis", ncores, quiet = quiet)
}

#' Minor axis length
#'
#' Computes the length of the minor axis (based on the bounding box) for each polygon.
#'
#' @inheritParams area
#' @return The sf object with an additional column `minor_axis`.
#' @references
#' Standard geometric measure.
#' @export
minor_axis <- function(sf_obj, ncores = 1, quiet = FALSE) {
  add_metric(sf_obj, calc_minor_axis, "minor_axis", ncores, quiet = quiet)
}

#' Perimeter of polygons
#'
#' Computes the perimeter of each polygon in the sf object.
#'
#' @inheritParams area
#' @return The sf object with an additional column `perimeter`.
#' @references
#' Standard geometric measure.
#' @export
perimeter <- function(sf_obj, ncores = 1, quiet = FALSE) {
  add_metric(sf_obj, calc_perimeter, "perimeter", ncores, quiet = quiet)
}

#' Perimeter-Area Ratio
#'
#' Computes the ratio of perimeter to area \eqn{P / A}. Lower values indicate
#' a more compact shape.
#'
#' @inheritParams area
#' @return The sf object with an additional column `perimeter_area_ratio`.
#' @references
#' McGarigal, K., Cushman, S. A., & Ene, E. (2012). FRAGSTATS v4.
#' @export
perimeter_area_ratio <- function(sf_obj, ncores = 1, quiet = FALSE) {
  add_metric(sf_obj, calc_perimeter_area_ratio, "perimeter_area_ratio", ncores, quiet = quiet)
}

#' Perimeter Index
#'
#' The perimeter index is the ratio of the perimeter of an equal-area circle
#' to the actual perimeter: \eqn{1 / \text{shape\_index}}. Values range from
#' 0 (very irregular) to 1 (perfect circle).
#'
#' @inheritParams area
#' @return The sf object with an additional column `perimeter_index`.
#' @references
#' McGarigal, K., Cushman, S. A., & Ene, E. (2012). FRAGSTATS v4.
#' @export
perimeter_index <- function(sf_obj, ncores = 1, quiet = FALSE) {
  add_metric(sf_obj, calc_perimeter_index, "perimeter_index", ncores, quiet = quiet)
}

#' Polsby‑Popper Compactness
#'
#' Computes the Polsby‑Popper index: \eqn{4\pi A / P^2}.
#' Values range from 0 (least compact) to 1 (perfect circle).
#'
#' @inheritParams area
#' @return The sf object with an additional column `polsby_popper`.
#' @references
#' Polsby, D. D., & Popper, R. D. (1991). The third criterion: Compactness as a
#' procedural safeguard against partisan gerrymandering.
#' *Yale Law & Policy Review*, 9(2), 301-353.
#' @export
polsby_popper <- function(sf_obj, ncores = 1, quiet = FALSE) {
  add_metric(sf_obj, calc_polsby_popper, "polsby_popper", ncores, quiet = quiet)
}

#' Proximity Index
#'
#' Ratio of the mean distance from the centroid in the equal‑area circle to
#' that in the polygon. Values close to 1 indicate that the polygon's internal
#' distances are similar to a circle.
#'
#' @inheritParams area
#' @param n_samples Number of random points to sample inside the polygon (default: 1000).
#' @return The sf object with an additional column `proximity`.
#' @references
#' McGarigal, K., Cushman, S. A., & Ene, E. (2012). FRAGSTATS v4.
#' @export
proximity <- function(sf_obj, ncores = 1, n_samples = 1000, quiet = FALSE) {
  previous_s2 <- disable_s2_if_geographic(sf_obj, "proximity")
  result <- add_metric(sf_obj, calc_proximity, "proximity", ncores, quiet = quiet, n_samples = n_samples)
  if (!is.null(previous_s2)) {
    sf_use_s2(previous_s2)
  }
  return(result)
}

#' Range Index
#'
#' Ratio of the diameter of the equal‑area circle to the diameter of the
#' minimum circumscribing circle. Values close to 1 indicate a shape that is
#' both compact and circular.
#'
#' @inheritParams area
#' @return The sf object with an additional column `range_index`.
#' @references
#' McGarigal, K., Cushman, S. A., & Ene, E. (2012). FRAGSTATS v4.
#' @export
range_index <- function(sf_obj, ncores = 1, quiet = FALSE) {
  add_metric(sf_obj, calc_range_index, "range_index", ncores, quiet = quiet)
}

#' Rectangularity (MABR)
#'
#' Computes the rectangularity index as the ratio of polygon area to the area
#' of its Minimum Area Bounding Rectangle (MABR). Values close to 1 indicate
#' a shape that is nearly rectangular.
#'
#' @inheritParams area
#' @return The sf object with an additional column `rectangularity`.
#' @references
#' Rosin, P. L. (2000). Measuring shape: ellipticity, rectangularity, and triangularity.
#' *Proceedings of the 15th International Conference on Pattern Recognition*.
#' Also used in FRAGSTATS.
#' @export
rectangularity <- function(sf_obj, ncores = 1, quiet = FALSE) {
  add_metric(sf_obj, calc_rectangularity, "rectangularity", ncores, quiet = quiet)
}

#' Reock Compactness
#'
#' Computes Reock index: area of polygon divided by area of its
#' minimum circumscribing circle.
#'
#' @inheritParams area
#' @references
#' Reock, E. C. (1961). A note: Measuring compactness as a requirement of
#' legislative apportionment. *Midwest Journal of Political Science*, 5(1), 70-74.
#' @return sf object with column `reock`.
#' @export
reock <- function(sf_obj, ncores = 1, quiet = FALSE) {
  add_metric(sf_obj, calc_reock, "reock", ncores, quiet = quiet)
}

#' Roughness
#'
#' Shape Index minus 1: \eqn{(P / (2 * \sqrt{\pi A})) - 1}.
#'
#' @inheritParams area
#' @return sf object with column `roughness`.
#' @references
#' McGarigal, K., Cushman, S. A., & Ene, E. (2012). FRAGSTATS v4.
#' @export
roughness <- function(sf_obj, ncores = 1, quiet = FALSE) {
  add_metric(sf_obj, calc_roughness, "roughness", ncores, quiet = quiet)
}

#' Schwartzberg Compactness
#'
#' Computes \eqn{1 / (P / (2\pi \sqrt{A/\pi}))}.
#'
#' @inheritParams area
#' @references
#' Schwartzberg, J. E. (1966). Reapportionment, gerrymanders, and the
#' 1964 elections. *The Annals of the American Academy of Political and Social Science*, 363(1), 54-68.
#' @return sf object with column `schwartzberg`.
#' @export
schwartzberg <- function(sf_obj, ncores = 1, quiet = FALSE) {
  add_metric(sf_obj, calc_schwartzberg, "schwartzberg", ncores, quiet = quiet)
}

#' Shape Index
#'
#' \eqn{P / (2 * \sqrt{\pi A})}. Minimum value is 1 for a circle.
#'
#' @inheritParams area
#' @return sf object with column `shape_index`.
#' @references
#' McGarigal, K., Cushman, S. A., & Ene, E. (2012). FRAGSTATS v4.
#' @export
shape_index <- function(sf_obj, ncores = 1, quiet = FALSE) {
  add_metric(sf_obj, calc_shape_index, "shape_index", ncores, quiet = quiet)
}

#' Skew Compactness
#'
#' Ratio of area of maximum inscribed circle to area of minimum circumscribing circle.
#'
#' @inheritParams area
#' @return sf object with column `skew`.
#' @references
#' McGarigal, K., Cushman, S. A., & Ene, E. (2012). FRAGSTATS v4.
#' @export
skew <- function(sf_obj, ncores = 1, quiet = FALSE) {
  previous_s2 <- disable_s2_if_geographic(sf_obj, "skew")
  result <- add_metric(sf_obj, calc_skew, "skew", ncores, quiet = quiet)
  if (!is.null(previous_s2)) {
    sf_use_s2(previous_s2)
  }
  return(result)
}

#' Sphericity
#'
#' Ratio of maximum inscribed circle radius to minimum circumscribing circle radius.
#'
#' @inheritParams area
#' @return sf object with column `sphericity`.
#' @references
#' McGarigal, K., Cushman, S. A., & Ene, E. (2012). FRAGSTATS v4.
#' @export
sphericity <- function(sf_obj, ncores = 1, quiet = FALSE) {
  previous_s2 <- disable_s2_if_geographic(sf_obj, "sphericity")
  result <- add_metric(sf_obj, calc_sphericity, "sphericity", ncores, quiet = quiet)
  if (!is.null(previous_s2)) {
    sf_use_s2(previous_s2)
  }
  return(result)
}

#' PGYRATIUS Index
#'
#' Computes the PGYRATIUS index: average distance of all vertices to the polygon centroid.
#'
#' @inheritParams area
#' @return The sf object with an additional column `p_gyradius`.
#' @references
#' Standard shape metric used in polygon analysis.
#' @export
p_gyradius <- function(sf_obj, ncores = 1, quiet = FALSE) {
  previous_s2 <- disable_s2_if_geographic(sf_obj, "p_gyradius")
  result <- add_metric(sf_obj, calc_p_gyradius, "p_gyradius", ncores, quiet = quiet)
  if (!is.null(previous_s2)) {
    sf_use_s2(previous_s2)
  }
  return(result)
}

#' POLRADIUS Index
#'
#' Computes the POLRADIUS index: maximum distance of a vertex to the polygon centroid.
#'
#' @inheritParams area
#' @return The sf object with an additional column `polradius`.
#' @references
#' Standard shape metric used in polygon analysis.
#' @export
polradius <- function(sf_obj, ncores = 1, quiet = FALSE) {
  previous_s2 <- disable_s2_if_geographic(sf_obj, "polradius")
  result <- add_metric(sf_obj, calc_polradius, "polradius", ncores, quiet = quiet)
  if (!is.null(previous_s2)) {
    sf_use_s2(previous_s2)
  }
  return(result)
}

#' PDENSITY Index
#'
#' Computes the PDENSITY index: \eqn{area / POLRADIUS}.
#'
#' @inheritParams area
#' @return The sf object with an additional column `p_density`.
#' @references
#' Standard shape metric used in polygon analysis.
#' @export
p_density <- function(sf_obj, ncores = 1, quiet = FALSE) {
  previous_s2 <- disable_s2_if_geographic(sf_obj, "p_density")
  result <- add_metric(sf_obj, calc_p_density, "p_density", ncores, quiet = quiet)
  if (!is.null(previous_s2)) {
    sf_use_s2(previous_s2)
  }
  return(result)
}

# =============================================================================
# SHARED BOUNDARY
# =============================================================================

#' Shared boundary length with a reference polygon
#'
#' Computes, for each polygon in \code{sf_obj}, the total length of its boundary
#' that is shared with a single reference polygon provided in \code{ref}.
#'
#' @param sf_obj An sf object with multiple polygon geometries.
#' @param ref An sf object containing exactly one polygon (the reference).
#' @param ncores Number of cores for parallel computation. Default is 1 (sequential).
#' @param quiet Logical: suppress messages about core usage.
#' @return The sf object with an additional column `shared_boundary`.
#' @references
#' Standard measure of spatial contiguity; used in landscape ecology and redistricting.
#' @export
shared_boundary <- function(sf_obj, ref, ncores = 1, quiet = FALSE) {
  check_sf(sf_obj)
  check_sf(ref, "ref")
  
  check_valid_geometry(sf_obj)
  check_valid_geometry(ref)
  
  check_polygon_geometry(sf_obj)
  check_polygon_geometry(ref, type = c("POLYGON", "MULTIPOLYGON"))
  
  previous_s2_main <- disable_s2_if_geographic(sf_obj, "shared_boundary")
  previous_s2_ref <- disable_s2_if_geographic(ref, "shared_boundary")
  
  if (nrow(ref) != 1) {
    stop("The reference object ('ref') must contain exactly one polygon.")
  }
  ref_poly <- ref[1, ]
  
  calc_shared <- function(poly, ref_poly) {
    inter <- st_intersection(poly, ref_poly)
    if (length(inter) == 0) return(0)
    geom_class <- class(inter)[1]
    if (grepl("LINESTRING|MULTILINESTRING", geom_class)) {
      return(as.numeric(st_length(inter)))
    } else if (grepl("GEOMETRYCOLLECTION", geom_class)) {
      lines <- st_collection_extract(inter, "LINESTRING")
      if (length(lines) > 0) {
        return(as.numeric(st_length(lines)))
      }
    }
    return(0)
  }
  
  add_metric_with_ref <- function(sf_obj, calc_func, ref_poly, col_name, ncores = 1, quiet = FALSE) {
    if (is.null(ncores) || ncores < 1) ncores <- 1
    ncores <- as.integer(ncores)
    max_cores <- future::availableCores()
    if (ncores > max_cores) {
      warning(sprintf("Requested %d cores, but only %d are available. Using %d cores.",
                      ncores, max_cores, max_cores - 1))
      ncores <- max(1, max_cores - 1)
    }
    
    if (ncores > 1 && !quiet) {
      message(sprintf("Using parallel processing with %d cores.", ncores))
    }
    
    if (ncores > 1) {
      future::plan("multisession", workers = ncores)
    }
    
    n_rows <- nrow(sf_obj)
    CHUNK_THRESHOLD <- 100000
    
    values <- if (n_rows > CHUNK_THRESHOLD) {
      chunk_size <- 5000
      n_chunks <- ceiling(n_rows / chunk_size)
      
      if (!quiet) {
        message(sprintf("Large dataset (%d polygons). Processing in %d chunks.", 
                        n_rows, n_chunks))
      }
      
      process_chunk <- function(chunk_idx) {
        start <- (chunk_idx - 1) * chunk_size + 1
        end <- min(chunk_idx * chunk_size, n_rows)
        chunk <- sf_obj[start:end, ]
        
        gc()
        
        chunk_results <- lapply(seq_len(nrow(chunk)), function(j) {
          tryCatch({
            calc_func(chunk[j, ], ref_poly)
          }, error = function(e) {
            warning(sprintf("Error calculating shared boundary for polygon %d: %s", 
                            start + j - 1, e$message))
            return(NA)
          })
        })
        
        # Liberar chunk da memória explicitamente
        rm(chunk)
        gc()
        
        return(chunk_results)
      }
      
      chunk_results_list <- tryCatch({
        future.apply::future_lapply(seq_len(n_chunks), process_chunk, 
                                    future.seed = TRUE, future.packages = "shapeMetrics")
      }, interrupt = function(e) {
        message("Processing interrupted by user.")
        stop(e)
      })
      
      unlist(chunk_results_list, recursive = FALSE)
      
    } else {
      # Modo normal: iterar por índices (sem split)
      tryCatch({
        future.apply::future_lapply(seq_len(n_rows), function(i) {
          tryCatch({
            calc_func(sf_obj[i, ], ref_poly)
          }, error = function(e) {
            warning(sprintf("Error calculating shared boundary for polygon %d: %s", i, e$message))
            return(NA)
          })
        }, future.seed = TRUE, future.packages = "shapeMetrics")
      }, interrupt = function(e) {
        message("Processing interrupted by user.")
        stop(e)
      })
    }
    
    gc()
    
    sf_obj[[col_name]] <- unlist(values)
    return(sf_obj)
  }
  
  result <- add_metric_with_ref(sf_obj, calc_shared, ref_poly, "shared_boundary", ncores, quiet = quiet)
  
  # Restore s2 (if it was altered)
  if (!is.null(previous_s2_main)) {
    sf_use_s2(previous_s2_main)
  }
  if (!is.null(previous_s2_ref)) {
    sf_use_s2(previous_s2_ref)
  }
  
  return(result)
}

# =============================================================================
# GROUP FUNCTIONS
# =============================================================================

#' Calculate Multiple Metrics
#'
#' Applies a list of metric functions to the sf object, adding each as a new column.
#' By default, includes all implemented shape/compactness metrics and geometric attributes.
#'
#' @param sf_obj An sf object.
#' @param metrics A named list of functions, a character vector of metric names, or \code{NULL}.
#'   If \code{NULL}, uses all available metrics. If a character vector, each element must
#'   correspond to a valid metric name (e.g., \code{c("area", "polsby_popper")}).
#' @param ncores Number of cores for each metric. Default is 1 (sequential).
#' @param progress Logical. If \code{TRUE}, displays a progress bar while processing metrics.
#'   Default is \code{FALSE}.
#' @param ... Additional arguments passed to each metric function.
#' @return The sf object with all metric columns added.
#' @references
#' The metrics are based on standard landscape ecology and spatial analysis literature.
#' @importFrom utils flush.console
#' @export
calc_multiple_metrics <- function(sf_obj, metrics = NULL, ncores = 1, progress = FALSE, ...) {
  check_sf(sf_obj)
  check_valid_geometry(sf_obj)
  check_polygon_geometry(sf_obj)
  
  if (ncores > 1) {
    message(sprintf("Using parallel processing with %d cores for all metrics.", ncores))
  }
  
  # --- Default metrics list ---
  default_metrics <- list(
    area                  = area,
    bounding_rect_area    = bounding_rect_area,
    box_reock             = box_reock,
    boyce_clark           = boyce_clark,
    circularity           = circularity,
    cohesion              = cohesion,
    compactness           = compactness,
    convex_hull_area      = convex_hull_area,
    convex_hull_compact   = convex_hull_compact,
    convex_hull_perimeter = convex_hull_perimeter,
    convexity             = convexity,
    detour                = detour,
    elongation            = elongation,
    equiv_rectangular     = equiv_rectangular,
    exchange              = exchange,
    fractal_dimension     = fractal_dimension,
    girth                 = girth,
    inscribed_radius      = inscribed_radius,
    major_axis            = major_axis,
    minor_axis            = minor_axis,
    p_density             = p_density,
    p_gyradius            = p_gyradius,
    perimeter             = perimeter,
    perimeter_area_ratio  = perimeter_area_ratio,
    perimeter_index       = perimeter_index,
    polradius             = polradius,
    polsby_popper         = polsby_popper,
    proximity             = proximity,
    range_index           = range_index,
    rectangularity        = rectangularity,
    reock                 = reock,
    roughness             = roughness,
    schwartzberg          = schwartzberg,
    shape_index           = shape_index,
    skew                  = skew,
    sphericity            = sphericity
  )
  valid_metric_names <- names(default_metrics)
  
  # --- Metrics validation ---
  if (is.null(metrics)) {
    metrics <- default_metrics
  } else {
    if (is.character(metrics)) {
      metric_names <- metrics
      valid_idx <- metric_names %in% valid_metric_names
      invalid_names <- metric_names[!valid_idx]
      
      if (length(invalid_names) > 0) {
        if (length(metric_names) == length(invalid_names)) {
          stop(sprintf(
            "None of the provided metric names are valid. Valid metrics are: %s",
            paste(valid_metric_names, collapse = ", ")
          ))
        } else {
          warning(sprintf(
            "The following metric(s) are not implemented and will be ignored: %s",
            paste(invalid_names, collapse = ", ")
          ))
          metric_names <- metric_names[valid_idx]
        }
      }
      if (length(metric_names) == 0) {
        stop("No valid metrics remaining after filtering.")
      }
      pkg_env <- asNamespace("shapeMetrics")
      metrics <- lapply(metric_names, function(name) {
        get(name, envir = pkg_env, mode = "function")
      })
      names(metrics) <- metric_names
    } else {
      if (is.null(names(metrics)) || any(names(metrics) == "")) {
        stop("The 'metrics' list must be a named list with metric names as names, or a character vector of metric names.")
      }
      invalid_names <- setdiff(names(metrics), valid_metric_names)
      if (length(invalid_names) > 0) {
        if (length(metrics) == length(invalid_names)) {
          stop(sprintf(
            "None of the provided metric names are valid. Valid metrics are: %s",
            paste(valid_metric_names, collapse = ", ")
          ))
        } else {
          warning(sprintf(
            "The following metric(s) are not implemented and will be ignored: %s",
            paste(invalid_names, collapse = ", ")
          ))
          metrics <- metrics[names(metrics) %in% valid_metric_names]
        }
      }
      if (length(metrics) == 0) {
        stop("No valid metrics remaining after filtering.")
      }
    }
    
    # shared_boundary handling
    if ("shared_boundary" %in% names(metrics)) {
      if (length(metrics) == 1) {
        stop("'shared_boundary' is not supported by calc_multiple_metrics because it requires an extra 'ref' argument. Use shared_boundary() directly.")
      } else {
        metrics <- metrics[names(metrics) != "shared_boundary"]
        warning("'shared_boundary' was removed from the metrics list because it requires an extra 'ref' argument.")
      }
    }
    if (length(metrics) == 0) {
      stop("No valid metrics remaining after removing 'shared_boundary'.")
    }
  }
  
  # --- Loop with progress (via cat) ---
  n_metrics <- length(metrics)
  metric_names <- names(metrics)
  
  if (progress && n_metrics > 0) {
    cat("\nProcessing metrics: 0%")
    flush.console()
  }
  
  for (i in seq_along(metrics)) {
    name <- metric_names[i]
    
    if (progress && n_metrics > 0) {
      pct <- round((i / n_metrics) * 100)
      cat(sprintf("\rProcessing metrics: %d%%", pct))
      flush.console()
    }
    
    sf_obj <- metrics[[name]](sf_obj, ncores = ncores, quiet = TRUE, ...)
  }
  
  if (progress && n_metrics > 0) {
    cat("\n")
  }
  
  return(sf_obj)
}

#' Compactness Metrics (redistmetrics style)
#'
#' Computes a subset of metrics focused on compactness as defined in the
#' \code{redistmetrics} package.
#'
#' @inheritParams calc_multiple_metrics
#' @return The sf object with compactness columns added.
#' @references
#' Based on the metrics described in the redistmetrics package.
#' @export
compactness_metrics <- function(sf_obj, ncores = 1, ...) {
  metrics <- list(
    polsby_popper     = polsby_popper,
    schwartzberg      = schwartzberg,
    reock             = reock,
    box_reock         = box_reock,
    convex_hull_compact = convex_hull_compact,
    skew              = skew,
    elongation        = elongation,
    boyce_clark       = boyce_clark,
    compactness       = compactness,
    perimeter_index   = perimeter_index,
    detour            = detour,
    girth             = girth,
    range_index       = range_index
  )
  calc_multiple_metrics(sf_obj, metrics, ncores, ...)
}

#' Shape Metrics (vectormetrics patch style)
#'
#' Computes a subset of metrics focused on shape, inspired by the patch metrics
#' in the \code{vectormetrics} package.
#'
#' @inheritParams calc_multiple_metrics
#' @return The sf object with shape metrics columns added.
#' @references
#' Inspired by the patch metrics described in the vectormetrics package.
#' @export
shape_metrics <- function(sf_obj, ncores = 1, ...) {
  metrics <- list(
    polsby_popper     = polsby_popper,
    circularity       = circularity,
    reock             = reock,
    box_reock         = box_reock,
    convex_hull_compact = convex_hull_compact,
    convexity         = convexity,
    sphericity        = sphericity,
    elongation        = elongation,
    fractal_dimension = fractal_dimension,
    shape_index       = shape_index,
    roughness         = roughness,
    compactness       = compactness,
    rectangularity    = rectangularity,
    perimeter_area_ratio = perimeter_area_ratio,
    perimeter_index   = perimeter_index,
    cohesion          = cohesion,
    detour            = detour,
    equiv_rectangular = equiv_rectangular,
    exchange          = exchange,
    girth             = girth,
    proximity         = proximity,
    range_index       = range_index
  )
  calc_multiple_metrics(sf_obj, metrics, ncores, ...)
}