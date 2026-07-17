# ====================================================================
# Tree trunk annulus erosion + Beta(u, K) profile integration
# Author: Michael W. Burnett, UC Santa Barbara
#
# INPUT:
#   - A GPKG of trunk cross-sections (one feature per tree) with columns:
#       tree_id (chr/int), geometry (MULTIPOLYGON/POLYGON),
#       perimeter_cm, hull_perimeter_cm, convexity,
#       area_cm2, hull_area_cm2, area_ratio, DBH_cm,
#       perimeter_to_area, deepest_concave_depth_cm, concave_depth_over_dbh,
#       [sapwood_depth_cm]  <-- required if mode = "sapwood"
#
# USER CONFIG:
#   - u, K for the Beta profile on z in (0,1)
#   - mode: "sapwood" if the beta distribution should be applied from d=0 to d=d_sw, or "full_radius" for d=0 to d=R_circ
#   - dx_mm erosion step (annulus thickness), 0.5 mm by default
#   - Amplitude scalar A, if real sapflow units are desired
#
# OUTPUT:
#   - One GPKG per tree_id in out_dir, layer name "layers", with attributes:
#       tree_id, layer (1..N), layer.distance.from.edge (mm),
#       layer.area.cm2, layer.area.mm2, layer.perimeter.cm,
#       swd.cm (if sapwood; NA otherwise),
#       relative.depth.to.swd (if sapwood) OR relative.depth.to.radius (if full),
#       relative.sfd  (unit-peak Beta(u,K) value at annulus mid-depth)
#
# Notes:
#   - "sapwood" mode requires the input GPKG has an additional column, sapwood_depth_cm, with sapwood depths for each tree!
#   - Parallelized via doParallel/foreach.
#
# SCRIPT 2 OF 4
#
# ====================================================================
#install.packages(c("sf","lwgeom","units","dplyr","tidyr","purrr","ggplot2","doParallel","foreach"))
suppressPackageStartupMessages({
  library(sf)
  library(lwgeom)
  library(units)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(ggplot2)
  library(doParallel)
  library(foreach)
})

# ----------------------------
# USER CONFIG
# ----------------------------
in_gpkg    <- "C:/trunks.gpkg" # change me
in_layer   <- sf::st_layers(in_gpkg)$name[1]
out_dir    <- "C:/eroded/" # change me

u          <- 0.348      # Beta distribution weighted mean (0,1)
K          <- 3.66       # Beta distribution concentration > 0
mode       <- "sapwood"  # "sapwood" or "full_radius"
dx_mm      <- 0.5        # annulus thickness in mm
crs_target <- 32601      # EPSG in meters
A          <- 1          # amplitude scale (1 for beta distribution)

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Validate u, K early
if (!(is.finite(u) && u > 0 && u < 1)) stop("u must be in (0,1).")
if (!(is.finite(K) && K > 0))          stop("K must be > 0.")

# ----------------------------
# Helpers
# ----------------------------

# Normalize to sfc MULTIPOLYGON
to_multipolygon <- function(g, crs = NA_crs_) {
  if (inherits(g, "sf"))  g <- st_geometry(g)
  if (inherits(g, "sfg")) g <- st_sfc(g, crs = crs)
  g <- st_make_valid(g)
  suppressWarnings(st_cast(g, "MULTIPOLYGON", warn = FALSE))
}

# Build annulus i (0-based): buffer(-i*dx) \ buffer(-(i+1)*dx)
one_ring <- function(g0, i, dx_m, crs) {
  d_out <- -i       * dx_m
  d_in  <- -(i + 1) * dx_m
  g_out <- suppressWarnings(st_buffer(g0, d_out))
  if (st_is_empty(g_out)) return(list(ok = FALSE))
  g_in  <- suppressWarnings(st_buffer(g0, d_in))
  ring  <- suppressWarnings(st_difference(g_out, g_in))
  
  g_out <- to_multipolygon(g_out, crs)
  g_in  <- to_multipolygon(g_in,  crs)
  ring  <- to_multipolygon(ring,  crs)
  if (st_is_empty(ring)) return(list(ok = FALSE))
  
  ring_u <- suppressWarnings(st_union(ring))
  ring_u <- to_multipolygon(ring_u, crs)
  list(ok = TRUE, ring = ring_u, g_out = g_out, g_in = g_in)
}

# Beta-pdf on z ∈ (0,1) parameterized by (u, K)
beta_pdf_muK <- function(z, u, K) {
  a <- u * K
  b <- (1 - u) * K
  z0 <- pmin(pmax(z, 1e-9), 1 - 1e-9)
  dbeta(z0, a, b)
}

# ----------------------------
# Read input & validate
# ----------------------------
polys <- st_read(in_gpkg, layer = in_layer, quiet = TRUE) |>
  st_transform(crs_target)

required_cols <- c(
  "tree_id","perimeter_cm","hull_perimeter_cm","convexity",
  "area_cm2","hull_area_cm2","area_ratio","DBH_cm",
  "perimeter_to_area","deepest_concave_depth_cm","concave_depth_over_dbh"
)
missing_cols <- setdiff(required_cols, names(polys))
if (length(missing_cols)) {
  stop("Missing required column(s): ", paste(missing_cols, collapse = ", "))
}
if (mode == "sapwood" && !"sapwood_depth_cm" %in% names(polys)) {
  stop("mode='sapwood' requires column 'sapwood_depth_cm' in the input GPKG.")
}

dx_m <- dx_mm / 1000   # mm -> m
crs  <- st_crs(polys)

# ----------------------------
# Parallel processing
# ----------------------------
cl <- makeCluster(max(1, parallel::detectCores() - 1))
on.exit(stopCluster(cl), add = TRUE)
registerDoParallel(cl)

res <- foreach(
  irow = seq_len(nrow(polys)),
  .export   = c("to_multipolygon","one_ring","beta_pdf_muK","u","K","mode","dx_m","dx_mm","A","out_dir","crs"),
  .packages = c("sf","lwgeom","units","dplyr","rlang")
) %dopar% {
  
  rec <- polys[irow, , drop = FALSE]
  g0  <- to_multipolygon(st_geometry(rec), crs = crs)
  if (length(g0) != 1) g0 <- g0[1]
  if (st_is_empty(g0)) return(NULL)
  
  # Depth normalizer for z ∈ (0,1)
  if (identical(mode, "sapwood")) {
    swd.cm <- as.numeric(rec$sapwood_depth_cm)
    if (!is.finite(swd.cm) || swd.cm <= 0) return(NULL)
    norm_cm <- swd.cm
    rel_depth_col <- "relative.depth.to.swd"
  } else if (identical(mode, "full_radius")) {
    radius_cm <- as.numeric(rec$DBH_cm) / 2
    if (!is.finite(radius_cm) || radius_cm <= 0) return(NULL)
    swd.cm <- NA_real_
    norm_cm <- radius_cm
    rel_depth_col <- "relative.depth.to.radius"
  } else {
    stop("mode must be either 'sapwood' or 'full_radius'.")
  }
  
  layers <- vector("list", 64L)  # pre-alloc a bit; will trim
  i <- 0L
  k <- 0L
  repeat {
    r <- one_ring(g0, i, dx_m = dx_m, crs = crs)
    if (!r$ok) break
    
    # Geometry-derived stats
    area_cm2 <- as.numeric(sum(st_area(r$g_out)) - sum(st_area(r$g_in))) * 1e4
    if (!is.finite(area_cm2) || area_cm2 <= 0) { i <- i + 1L; next }
    area_cm2_hp <- round(area_cm2, 6)
    area_mm2    <- round(area_cm2 * 100, 0)     # 1 cm^2 = 100 mm^2
    perim_cm    <- as.numeric(lwgeom::st_perimeter(r$g_out)) * 100
    
    layer_idx      <- i + 1L
    layer_dist_mm  <- i * dx_mm
    mid_depth_cm   <- (layer_dist_mm + dx_mm/2) / 10
    z              <- mid_depth_cm / norm_cm
    z              <- max(min(z, 1 - 1e-9), 1e-9)
    
    # Raw pdf scaled by amplitude A
    rel_sfd <- A * beta_pdf_muK(z, u = u, K = K)
    
    # Relative sapflow for this annulus (use unrounded area)
    layer_relative_sapflow <- rel_sfd * area_cm2
    
    # Attribute row (with dynamic relative-depth column)
    dat_row <- rec |>
      st_drop_geometry() |>
      mutate(
        layer = layer_idx,
        layer.distance.from.edge = layer_dist_mm,
        layer.area.cm2 = area_cm2_hp,
        layer.area.mm2 = area_mm2,
        layer.perimeter.cm = perim_cm,
        swd.cm = if (identical(mode, "sapwood")) swd.cm else NA_real_,
        relative.sfd = rel_sfd,
        layer.relative.sapflow = layer_relative_sapflow
      )
    dat_row[[rel_depth_col]] <- z
    
    k <- k + 1L
    layers[[k]] <- st_sf(dat_row, geometry = st_geometry(r$ring))
    i <- i + 1L
  }
  
  if (k == 0L) return(NULL)
  layers <- layers[seq_len(k)]
  out <- do.call(rbind, layers)
  st_crs(out) <- crs
  
  out_path <- file.path(out_dir, sprintf("%s_eroded.gpkg", out$tree_id[1]))
  if (file.exists(out_path)) unlink(out_path)
  st_write(out, out_path, layer = "layers", quiet = TRUE)
  
  out[1:min(3, nrow(out)), c("tree_id","layer","layer.area.cm2","layer.perimeter.cm")]
}

# Optional: print summaries
summaries <- purrr::compact(res)
if (length(summaries)) {
  message("Wrote ", length(summaries), " trunks to: ", out_dir)
  print(summaries[[1]])
} else {
  message("No outputs produced (check inputs and mode/columns).")
}
