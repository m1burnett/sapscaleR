# =============================================================================
# Compute geometric attributes for tree trunk cross-sections (GPKG)
# - Reads a spatial layer with polygons and a "tree_id" field
# - Ensures valid geometries and metric CRS
# - Computes perimeter, convex hull perimeter, area, convex hull area, DBH, convexity, and deepest concavity depth
# - Writes a new GPKG with appended attributes
# - OPTIONAL: Uncomment the two code blocks for plotting to verify that the concavity depth measurements are correct
# - OPTIONAL: Uncomment the "sapwood_depth_path" user input line and "sapwood depth (cm)" section if you have sapwood depth data for each tree and later want to apply sapflux density profiles across active sapwood depths only (uncommented by default).
#
# SCRIPT 1 OF 4
#
# Author: Michael W. Burnett, UC Santa Barbara
# =============================================================================

# ---- Dependencies ------------------------------------------------------------
#install.packages(c("sf","lwgeom","dplyr","units"))
library(sf)
library(lwgeom)
library(dplyr)
library(units)

# ---- User inputs -------------------------------------------------------------
gpkg_path  <- "C:/smoothed_polygons_aligned_clean.gpkg" # <- change me
layer_name <- NULL            # set to a string if not the first layer
out_path   <- "C:/trunks.gpkg" # <- change me
crs_target <- 32601           # UTM zone 1N (meters). Change if needed.
sapwood_depth_path <- "C:/swd.csv"  # <- change me

# ---- Helpers -----------------------------------------------------------------

# Safe perimeter function across lwgeom versions
.perimeter_m <- function(g) {
  if ("st_perimeter" %in% getNamespaceExports("lwgeom")) {
    lwgeom::st_perimeter(g)
  } else {
    lwgeom::st_perimeter_lwgeom(g)  # fallback for older lwgeom
  }
}

# Return (x,y) of first coordinate of an sf point/geometry
.xy_first <- function(g) {
  cxy <- sf::st_coordinates(g)
  if (is.null(dim(cxy))) as.numeric(cxy[1:2]) else as.numeric(cxy[1, 1:2])
}

# Compute deepest concavity depth (cm) by casting inward rays from convex hull
deepest_concave_depth_cm_one <- function(poly,
                                         step_m = 0.005,   # sampling step along hull (m)
                                         tree_id = NA_character_,
                                         do_plot = FALSE) {
  stopifnot(inherits(poly, "sfc") || inherits(poly, "sfg"))
  crs_here <- sf::st_crs(poly)
  
  # Clean topology
  poly <- sf::st_buffer(poly, 0)
  if (sf::st_is_empty(poly)) return(0)
  
  # Convex hull
  hull <- sf::st_convex_hull(poly)
  
  # If convex (areas equal), depth = 0
  if (isTRUE(all.equal(as.numeric(sf::st_area(hull)),
                       as.numeric(sf::st_area(poly))))) {
    return(0)
  }
  
  # Sample points along hull boundary
  hull_line <- sf::st_boundary(hull)
  L <- as.numeric(sf::st_length(hull_line))
  if (!is.finite(L) || L == 0) return(0)
  n  <- max(4, ceiling(L / step_m) + 1)  # ensure wrap-around coverage
  pts <- sf::st_line_sample(hull_line, n = n, type = "regular") |> sf::st_cast("POINT")
  C   <- sf::st_coordinates(pts)[, 1:2, drop = TRUE]
  N   <- nrow(C)
  
  # Ray length: 2 × bbox diagonal of hull
  bb <- sf::st_bbox(hull)
  ray_len <- 2 * sqrt((bb$xmax - bb$xmin)^2 + (bb$ymax - bb$ymin)^2)
  
  # Centroid of hull for "inward" decision
  cen <- .xy_first(sf::st_centroid(hull))
  
  # Polygon boundary for intersections
  poly_boundary <- sf::st_boundary(poly)
  
  # Trackers
  max_depth_m <- 0
  p_start_max <- NULL
  ray_max     <- NULL
  interP_max  <- NULL
  
  # Optional base plot (un-comment this if you want to verify that the concavity measurements are correct)
  # if (do_plot) {
  #   plot(sf::st_geometry(poly), col = adjustcolor("grey70", 0.4), border = "grey30", asp = 1)
  #   plot(sf::st_geometry(hull), add = TRUE, border = "red", lwd = 2)
  # }
  
  for (i in seq_len(N)) {
    # Neighbor-based tangent (wrap-around)
    i_prev <- if (i == 1) N else i - 1
    i_next <- if (i == N) 1 else i + 1
    c_prev <- C[i_prev, ]; c_next <- C[i_next, ]; c0 <- C[i, ]
    
    # Tangent and normals
    tv <- c_next - c_prev
    tv_len <- sqrt(sum(tv^2))
    if (!is.finite(tv_len) || tv_len == 0) next
    t_hat <- tv / tv_len
    n_left  <- c(-t_hat[2],  t_hat[1])
    n_right <- c( t_hat[2], -t_hat[1])
    
    # Inward = side whose normal points more toward hull centroid
    v_cent <- cen - c0
    n_in <- if (sum(n_left * v_cent) >= sum(n_right * v_cent)) n_left else n_right
    
    # Start point: nudge slightly away from centroid to avoid grazing 0-distance hits
    p_start <- sf::st_sfc(sf::st_point(c0), crs = crs_here)
    cent_xy <- cen
    p0      <- .xy_first(p_start)
    v       <- p0 - cent_xy
    v       <- v / sqrt(sum(v^2))
    eps     <- 0.003                     # 3 mm nudge; negligible vs. rays
    p_start_nudged <- sf::st_sfc(sf::st_point(p0 + v * eps), crs = crs_here)
    
    p_end <- sf::st_sfc(sf::st_point(c0 + n_in * ray_len), crs = crs_here)
    ray   <- sf::st_sfc(sf::st_linestring(rbind(.xy_first(p_start_nudged), .xy_first(p_end))), crs = crs_here)
    
    # Intersections of ray with polygon boundary
    interP <- suppressWarnings(sf::st_intersection(ray, poly_boundary))
    if (sf::st_is_empty(interP)) next
    pts_hit <- suppressWarnings(sf::st_collection_extract(interP, "POINT"))
    if (length(pts_hit) == 0) next
    
    # Distances from p_start to each hit (along the ray)
    cc_hit <- sf::st_coordinates(pts_hit)[, 1:2, drop = FALSE]
    p0     <- .xy_first(p_start_nudged)
    dists  <- sqrt((cc_hit[,1] - p0[1])^2 + (cc_hit[,2] - p0[2])^2)
    if (!length(dists)) next
    
    # Local depth = minimum positive distance (includes near-zero for convex contact)
    d_min <- min(dists, na.rm = TRUE)
    if (is.finite(d_min) && d_min > max_depth_m) {
      max_depth_m <- d_min
      p_start_max <- p_start_nudged
      ray_max     <- ray
      interP_max  <- pts_hit
    }
  }
  
  # Optional overlay of chosen ray & points (un-comment this if you want to verify that the concavity measurements are correct)
  # if (do_plot && !is.null(ray_max)) {
  #   points(sf::st_geometry(p_start_max), pch = 21, bg = "gold", col = "black", cex = 1.2)
  #   plot(sf::st_geometry(ray_max), add = TRUE, col = "dodgerblue", lwd = 2)
  #   plot(sf::st_geometry(interP_max), add = TRUE, pch = 21, bg = "red", col = "black", cex = 1.1)
  #   ttl <- sprintf("Tree %s — Deepest concavity = %.1f cm",
  #                  ifelse(is.na(tree_id), "<?>", tree_id), max_depth_m * 100)
  #   title(ttl, font.main = 2)
  # }
  
  # Return centimeters
  max_depth_m * 100
}


# ---- Read layer --------------------------------------------------------------
if (is.null(layer_name)) {
  layer_name <- st_layers(gpkg_path)$name[1]
}
x <- st_read(gpkg_path, layer = layer_name, quiet = TRUE)

# ---- Keep only tree_id + geometry, ensure validity --------------------------------
stopifnot("tree_id" %in% names(x))
x <- st_make_valid(x)
x <- st_sf(tree_id = x$tree_id, geometry = st_geometry(x))

# ---- Ensure metric CRS (transform if input is lon/lat or not in target CRS) ---
crs_in <- st_crs(x)
if (is.na(crs_in)) {
  message("Input has no CRS; assigning target CRS EPSG:", crs_target)
  st_crs(x) <- crs_target
} else if (isTRUE(st_is_longlat(x)) || (!is.na(crs_in$epsg) && crs_in$epsg != crs_target)) {
  message("Transforming to EPSG:", crs_target, " for accurate metric measurements.")
  x <- st_transform(x, crs_target)
}

# ---- Geometric attributes (meters) -------------------------------------------
perim_m      <- .perimeter_m(x)               # perimeter (includes holes if present)
hull_geom    <- st_convex_hull(x)
hull_perim_m <- .perimeter_m(hull_geom)
area_m2      <- st_area(x)
hull_area_m2 <- st_area(hull_geom)

# ---- Convert to centimeters and append ---------------------------------------
x$perimeter_cm       <- as.numeric(set_units(perim_m, "cm"))
x$hull_perimeter_cm  <- as.numeric(set_units(hull_perim_m, "cm"))
x$convexity          <- x$hull_perimeter_cm / x$perimeter_cm   # (0,1], 1 = perfectly convex/circular
x$area_cm2           <- as.numeric(set_units(area_m2, "cm^2"))
x$hull_area_cm2      <- as.numeric(set_units(hull_area_m2, "cm^2"))
x$area_ratio         <- x$area_cm2 / x$hull_area_cm2            # ≥ 1 if shape is concave wrt hull area
x$DBH_cm             <- x$hull_perimeter_cm / pi                # equivalent-circle diameter
x$perimeter_to_area  <- x$perimeter_cm / x$area_cm2

# ---- Deepest concavity depth (cm) --------------------------------------------
# Vectorize across features — IMPORTANT: pass sfc geometry with st_geometry(x)[i]
x$deepest_concave_depth_cm <- vapply(
  seq_len(nrow(x)),
  function(i) deepest_concave_depth_cm_one(
    sf::st_geometry(x)[i],
    step_m  = 0.005,
    tree_id = x$tree_id[i],
    do_plot = do_plot
  ),
  numeric(1)
)

# Derived scalars
x$concave_depth_over_dbh       <- x$deepest_concave_depth_cm / x$DBH_cm

# ---- sapwood depth (cm) - UNCOMMENT TO ADD SWD IF YOUR RADIAL SAPFLUX PROFILES SPAN ACTIVE SAPWOOD ONLY --------------------------------------------
# inserts sapwood depth data from other file
swd <- read.csv(sapwood_depth_path)
swd$tree_id <- as.integer(swd$tree_id)
x <- left_join(x, swd[,c(2,3)], by="tree_id")

# relative sapwood depth r_sapwood
x <- x %>%
  mutate(
    sap_cm = suppressWarnings(as.numeric(sapwood_depth_cm)),
    rad_cm = suppressWarnings(as.numeric(DBH_cm)) / 2,
    # relative sapwood depth; if depth > radius, use 1
    r_sapwood = ifelse(is.finite(sap_cm) & is.finite(rad_cm) & rad_cm > 0,
                     pmin(sap_cm / rad_cm, 1),    # clamp at 1 = deeper than radius
                     NA_real_)
  ) %>%
  select(-sap_cm, -rad_cm)

# ---- Write out ---------------------------------------------------------------
st_write(x, out_path, layer = layer_name, delete_dsn = TRUE, quiet = TRUE)
message("Wrote: ", out_path)
