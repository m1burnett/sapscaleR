# =============================================================================
# Compute r_steiner and r_in for each trunk polygon
#
# SCRIPT 3 OF 4
#
# Author: Michael W. Burnett, UC Santa Barbara
# =============================================================================

#install.packages(c("sf","dplyr","tidyr","purrr","tibble"))
suppressPackageStartupMessages({
  library(sf)
  library(dplyr)
  library(purrr)
  library(tidyr)
  library(tibble)
})


# --- Steiner depth (r_steiner) for each polygon, as estimated by an absolute deviation from Steiner-predicted relative annulus area of 0.02 ----------------

gpkg_path <- "C:/trunks.gpkg" # change me
eroded_dir <- "C:/eroded/" # change me


# ---- helpers ----

# Build g(r) from inward buffers.
# r = t / R_h, where R_h = P_hull / (2*pi).  g(t) = (-dA/dt) / P_hull.
g_curve_from_poly <- function(geom, R_h_cm, P_hull_cm,
                              step_cm = 0.25, rc_max = 1.0,
                              target_pts = 300, min_step_cm = 0.02) {
  geom <- tryCatch(sf::st_make_valid(geom), error = function(e) geom)
  geom <- sf::st_cast(geom, "MULTIPOLYGON", warn = FALSE)
  
  t_max_cm <- rc_max * R_h_cm
  step_eff <- max(min_step_cm, min(step_cm, t_max_cm / target_pts))
  ts_cm <- seq(0, t_max_cm, by = step_eff)
  if (length(ts_cm) < 3) return(NULL)
  
  # Areas for inward offsets
  A <- numeric(length(ts_cm)); ok <- rep(TRUE, length(ts_cm))
  for (i in seq_along(ts_cm)) {
    buf <- try(sf::st_buffer(geom, -ts_cm[i] / 100), silent = TRUE) # cm -> m
    if (inherits(buf, "try-error") || any(sf::st_is_empty(buf))) { ok[i] <- FALSE; break }
    A[i] <- as.numeric(sum(sf::st_area(buf))) # m^2
  }
  last_ok <- max(which(ok))
  if (!is.finite(last_ok) || last_ok < 3) return(NULL)
  
  ts_cm <- ts_cm[1:last_ok]; A <- A[1:last_ok]
  n <- length(ts_cm)
  
  # Central-diff derivative in meters
  dt_m    <- (ts_cm[3:n] - ts_cm[1:(n-2)]) / 100
  dA_dt_m <- - (A[3:n] - A[1:(n-2)]) / dt_m   # ≈ perimeter at offset, in meters
  t_mid_cm <- ts_cm[2:(n-1)]
  
  if (!is.finite(dA_dt_m[1]) || dA_dt_m[1] <= 0) return(NULL)
  
  P_hull_m <- P_hull_cm / 100
  g  <- pmax(dA_dt_m / P_hull_m, 0)         
  r  <- t_mid_cm / R_h_cm
  
  distinct(tibble(r_c = r, g = g), r_c, .keep_all = TRUE)
}

# Detect r_steiner via absolute vertical gap: first sustained r where (g_Steiner - g) >= abs_gap
detect_r_steiner <- function(curve, alpha,
                             abs_gap = 0.02, min_rc = 0.02,
                             ma_width_rc = 0.02, run_len_rc = 0.01) {
  if (is.null(curve) || nrow(curve) < 3) return(NA_real_)
  df <- dplyr::arrange(curve, r_c)
  
  # Regular grid for stable indexing
  r <- seq(max(min(df$r_c), 0), min(max(df$r_c), 1), length.out = 400)
  g <- approx(df$r_c, df$g, xout = r, ties = "ordered")$y
  
  # ---- Steiner reference: gS = PR - r = (1/alpha) - r ----
  PR  <- 1/alpha
  gS  <- pmax(0, PR - r)
  
  gap <- gS - g
  
  # Moving average over ~ma_width_rc in r
  dr <- stats::median(diff(r), na.rm = TRUE)
  k  <- max(3, round(ma_width_rc / dr))
  gap_ma <- as.numeric(stats::filter(gap, rep(1/k, k), sides = 1))
  
  # Find first sustained run >= abs_gap for ~run_len_rc
  L <- max(2, round(run_len_rc / dr))
  ok <- which(r >= min_rc & is.finite(gap_ma))
  first <- NA_integer_
  if (length(ok)) {
    above <- gap_ma >= abs_gap
    for (i in ok) {
      j2 <- i + L - 1L
      if (j2 <= length(above) && all(above[i:j2], na.rm = FALSE)) { first <- i; break }
    }
  }
  
  if (is.na(first)) {
    return(max(r, na.rm = TRUE))  # never exceeded cutoff
  } else {
    # Linear interpolation to approximate crossing
    j <- max(which(gap_ma[seq_len(first)] < abs_gap), na.rm = TRUE)
    if (is.finite(j) && j >= 1L) {
      x0 <- r[j];     y0 <- gap_ma[j] - abs_gap
      x1 <- r[first]; y1 <- gap_ma[first] - abs_gap
      if (is.finite(y0) && is.finite(y1) && (y1 - y0) != 0)
        return(x0 + (x1 - x0) * (-y0) / (y1 - y0))
    }
    return(r[first])
  }
}


# ---- read gpkg ----
layer_name <- sf::st_layers(gpkg_path)$name[1]
dat <- sf::st_read(gpkg_path, layer = layer_name, quiet = TRUE)

# required columns check
need <- c("perimeter_cm", "hull_perimeter_cm")
stopifnot(all(need %in% names(dat)))

# required columns check
need <- c("perimeter_cm", "hull_perimeter_cm", "convexity")
stopifnot(all(need %in% names(dat)))

# pull alpha from stored convexity; fallback to hull_perimeter/perimeter if needed
alpha_from_col <- suppressWarnings(as.numeric(dat$convexity))
alpha_fallback <- dat$hull_perimeter_cm / dat$perimeter_cm
alpha_use <- ifelse(is.finite(alpha_from_col) & alpha_from_col > 0, alpha_from_col, alpha_fallback)

# precompute per-row constants
meta <- dat |>
  mutate(
    R_h_cm = hull_perimeter_cm / (2 * pi),   # radius of circle with circumference = hull perimeter
    alpha  = alpha_use,                     
    .row   = dplyr::row_number()
  ) |>
  st_drop_geometry()

# compute r_steiner per polygon
rvals <- pmap_dbl(
  list(idx = meta$.row,
       R_h = meta$R_h_cm,
       P_h = meta$hull_perimeter_cm,
       a   = meta$alpha),
  function(idx, R_h, P_h, a) {
    geom <- sf::st_geometry(dat)[idx]
    # extract polygons from geometry collections if any
    if (inherits(geom[[1]], "GEOMETRYCOLLECTION")) {
      geom <- sf::st_collection_extract(geom, "POLYGON")
    }
    curve <- tryCatch(
      g_curve_from_poly(geom, R_h_cm = R_h, P_hull_cm = P_h,
                        step_cm = 0.25, rc_max = 1.0, target_pts = 300, min_step_cm = 0.02),
      error = function(e) NULL
    )
    detect_r_steiner(curve, alpha = a,
                     abs_gap = 0.02, min_rc = 0.02,
                     ma_width_rc = 0.02, run_len_rc = 0.01)
  }
)

# attach and write back
dat$r_steiner <- rvals

# how much does sapwood depth exceed r_steiner?
# if these values are negative, active sapwood is all shallower than r_steiner.
# if positive, active sapwood exists beyond r_steiner so SF_steiner must be applied with caution.
if ("r_sapwood" %in% names(dat)) {
  dat <- dat %>% mutate(
      swd_minus_r_steiner = r_sapwood - r_steiner
    )
  message("Added 'swd_minus_r_steiner' = min(sapwood_depth_cm/(DBH_cm/2), 1) - r_steiner.")
} else {
  message("Column 'sapwood_depth_cm' not found; skipping 'swd_minus_r_steiner'.")
}

dat$incircle_R_cm <- NA_real_
dat$R_cm          <- as.numeric(dat$DBH_cm) / 2
dat$r_in          <- NA_real_

for (i in seq_len(nrow(dat))) {
  tree_id_chr <- trimws(as.character(dat$tree_id[i]))
  gpkg_path_eroded   <- file.path(eroded_dir, sprintf("%s_eroded.gpkg", tree_id_chr))
  
  if (!file.exists(gpkg_path_eroded)) {
    warning("Missing file for tree_id=", tree_id_chr, ": ", gpkg_path_eroded)
    next
  }
  
  # Pick layer (prefer "layers")
  lyr_names <- sf::st_layers(gpkg_path_eroded)$name
  lyr <- if ("layers" %in% lyr_names) "layers" else lyr_names[1]
  
  x <- suppressMessages(sf::st_read(gpkg_path_eroded, layer = lyr, quiet = TRUE))
  if (!"layer.distance.from.edge" %in% names(x)) {
    warning("No 'layer.distance.from.edge' in ", basename(gpkg_path_eroded))
    next
  }
  
  # Max distance (mm) → cm
  max_dist_mm <- suppressWarnings(max(as.numeric(x$layer.distance.from.edge), na.rm = TRUE))
  if (!is.finite(max_dist_mm)) next
  
  incircle_R_cm <- max_dist_mm / 10  # mm → cm
  dat$incircle_R_cm[i] <- incircle_R_cm
  
  R_cm_i <- dat$R_cm[i]
  dat$r_in[i] <- if (is.finite(R_cm_i) && R_cm_i > 0) incircle_R_cm / R_cm_i else NA_real_
}


sf::st_write(dat, dsn = gpkg_path, layer = layer_name, delete_layer = TRUE, quiet = FALSE)
cat("Done. Wrote r_steiner to layer:", layer_name, "\n")
