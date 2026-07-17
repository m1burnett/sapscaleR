# ====================================================================
# Whole-tree sapflow aggregation (from eroded annuli) + SF_circ, SF_steiner, and SF_corr estimation
#
# SCRIPT 4 OF 4
#
#	- By default, SF_corr is estimated not with the r_steiner and r_in values measured by the previous script, but instead with modeled values from Eqs. 10-13 in the article. Measured values can be used instead by changing "use_regression_for_r" to FALSE; regression coefficients can also be modified in the "toggle & coefficients" section. The script uses Eqs. 10-13 by default because in the article we were examining the overall uncertainty in SF_corr if the models are used for r_steiner and r_in, as they would be in a new tree. Eqs. 10-13 were fit using the data from the previous script.
# - Mode can be changed from default "sapwood" to "full_radius" if the beta sapflux distribution reflects the entire range of depths from cambium to R_circ rather than the range of depths from cambium to sapwood depth. Stay consistent with the mode used in trunk-erosion.R
# - Scaling parameter A should be the same in this script as in Script #2 (1 by default in both).
#
# Author: Michael W. Burnett, UC Santa Barbara
# ====================================================================

#install.packages(c("sf","dplyr","stringr","purrr","units"))
suppressPackageStartupMessages({
  library(sf)
  library(dplyr)
  library(purrr)
  library(stringr)
  library(units)
})

# ----------------------------
# USER CONFIG
# ----------------------------
eroded_dir <- "C:/eroded/" # change me

master_gpkg  <- "C:/trunks.gpkg" # change me
master_layer <- sf::st_layers(master_gpkg)$name[1]   # will overwrite this layer

mode <- "sapwood"      # "sapwood" or "full_radius"
u    <- 0.348          # beta distribution weighted mean in (0,1)
K    <- 3.66           # total shape mass -> alpha=u*K, beta=(1-u)*K
A    <- 1              # amplitude multiplier; A=1 => pdf

# === toggle & coefficients ================
use_regression_for_r <- TRUE   # <- turn off to instead use measured r_in and r_steiner
# r_steiner model: logit(mu) = b0_st + b1_st * logit(c_p)
b0_st <- -2.885; b1_st <-  0.634
# r_in model:     logit(mu) = b0_in + b1_in * logit(c_p) + b2_in * h_rel
b0_in <-  0.303; b1_in <-  0.297; b2_in <- -1.767
# above values are the functions fit to Pisonia grandis

# ----------------------------
# Validations
# ----------------------------
if (!(is.finite(u) && u > 0 && u < 1)) stop("u must be in (0,1).")
if (!(is.finite(K) && K > 0))          stop("K must be > 0.")

# ----------------------------
# Helpers
# ----------------------------

eps      <- 1e-9
invlogit <- plogis
logit    <- function(x) log(pmin(pmax(x, eps), 1-eps) / (1 - pmin(pmax(x, eps), 1-eps)))
clamp01  <- function(x) pmin(pmax(x, 0), 1)

# Try to derive h_rel (= DCD/DBH) from available columns, else NA (handled in mutate)
derive_h_rel <- function(row) {
  nms <- names(row)
  if ("h_rel" %in% nms && is.finite(row$h_rel)) return(as.numeric(row$h_rel))
  if ("concave_depth_over_dbh" %in% nms && is.finite(row$concave_depth_over_dbh))
    return(as.numeric(row$concave_depth_over_dbh))
  if (all(c("deepest_concave_depth_cm","DBH_cm") %in% nms)) {
    d <- as.numeric(row$deepest_concave_depth_cm)
    D <- as.numeric(row$DBH_cm)
    if (is.finite(d) && is.finite(D) && D > 0) return(d / D)
  }
  if (all(c("DCD_cm","DBH_cm") %in% nms)) {
    d <- as.numeric(row$DCD_cm); D <- as.numeric(row$DBH_cm)
    if (is.finite(d) && is.finite(D) && D > 0) return(d / D)
  }
  NA_real_
}

predict_r_steiner_mu <- function(cp) {
  invlogit(b0_st + b1_st * logit(cp))
}

predict_r_in_mu <- function(cp, hrel) {
  invlogit(b0_in + b1_in * logit(cp) + b2_in * hrel)
}

sum_tree_sapflow <- function(gpkg_path) {
  lyr_names <- sf::st_layers(gpkg_path)$name
  lyr <- if ("layers" %in% lyr_names) "layers" else lyr_names[1]
  x <- suppressMessages(sf::st_read(gpkg_path, layer = lyr, quiet = TRUE))
  tree_id <- if ("tree_id" %in% names(x)) as.character(unique(x$tree_id))[1]
  if (!"layer.relative.sapflow" %in% names(x)) {
    stop("Missing 'layer.relative.sapflow' in ", basename(gpkg_path))
  }
  sapflow_ref <- sum(as.numeric(x$layer.relative.sapflow), na.rm = TRUE)
  tibble::tibble(tree_id = as.character(tree_id), sapflow_ref = sapflow_ref)
}

sapflow_circle_closed_form <- function(R_cm, p, u, A) {
  p <- pmin(pmax(p, 0), 1)
  2 * pi * (R_cm^2) * p * A * (1 - u * p)
}

sf_corr_beta <- function(R_cm, p, cp, r_steiner, r_in, alpha, beta, A = 1) {
  if (!is.finite(R_cm) || !is.finite(p) || !is.finite(cp) ||
      !is.finite(r_steiner) || !is.finite(r_in) ||
      !is.finite(alpha) || !is.finite(beta) ||
      R_cm <= 0 || p <= 0 || cp <= 0 || alpha <= 0 || beta <= 0) {
    return(NA_real_)
  }
  
  # z coordinates
  z_s_raw  <- r_steiner / p
  z_in_raw <- r_in / p
  
  # integration bounds are clamped to active sapwood domain [0, 1]
  z_s     <- max(0, min(1, z_s_raw))
  z_upper <- max(0, min(1, z_in_raw))
  
  F  <- function(x) pbeta(x, alpha,     beta)
  G1 <- function(x) pbeta(x, alpha + 1, beta)
  mu <- alpha / (alpha + beta)
  
  # I1 = ∫_0^z_s (1 - cp*p*z) rho(z) dz
  I1 <- F(z_s) - cp * p * mu * G1(z_s)
  
  # I2 = ∫_z_s^z_upper (1 - cp*r_steiner) *
  #      (r_in - p*z)/(r_in - r_steiner) * rho(z) dz
  I2 <- 0
  if (z_upper > z_s) {
    denom <- max(r_in - r_steiner, .Machine$double.eps)
    fac   <- (1 - cp * r_steiner) / denom
    
    dF <- F(z_upper) - F(z_s)
    dG <- G1(z_upper) - G1(z_s)
    
    # ∫ (r_in - p*z) rho(z) dz
    I2 <- fac * (r_in * dF - p * mu * dG)
  }
  
  pref <- 2 * pi * R_cm^2 * p * A / cp
  as.numeric(pref * (I1 + I2))
}

# ----------------------------
# 1) Aggregate sapflow_ref from eroded per-tree GPKGs
# ----------------------------
gpkg_files <- list.files(eroded_dir, pattern = "\\.gpkg$", full.names = TRUE)
if (!length(gpkg_files)) stop("No GPKG files found in: ", eroded_dir)

sapflow_ref_df <- gpkg_files %>%
  map_dfr(sum_tree_sapflow) %>%
  distinct(tree_id, .keep_all = TRUE)

# ----------------------------
# 2) Read master trunks and join sapflow_ref
# ----------------------------
trunks <- suppressMessages(sf::st_read(master_gpkg, layer = master_layer, quiet = TRUE))
trunks <- trunks %>% mutate(tree_id = trimws(as.character(tree_id)))
sapflow_ref_df <- sapflow_ref_df %>% mutate(tree_id = trimws(as.character(tree_id)))
trunks <- trunks %>% left_join(sapflow_ref_df, by = "tree_id")

# If the join produced sapflow_ref.x / sapflow_ref.y, coalesce into sapflow_ref
if (!"sapflow_ref" %in% names(trunks)) {
  if (all(c("sapflow_ref.x","sapflow_ref.y") %in% names(trunks))) {
    trunks <- trunks %>%
      mutate(sapflow_ref = dplyr::coalesce(as.numeric(sapflow_ref.x),
                                           as.numeric(sapflow_ref.y))) %>%
      select(-sapflow_ref.x, -sapflow_ref.y)
  }
}

n_joined <- sum(!is.na(trunks$sapflow_ref))
n_total  <- nrow(trunks)
if (n_joined == 0) {
  warning("No sapflow_ref values joined. Check eroded_dir, file names, and tree_id.")
} else if (n_joined < n_total) {
  warning(sprintf("sapflow_ref joined for %d/%d trees. Some tree_id keys did not match.",
                  n_joined, n_total))
}

if (!"sapflow_ref" %in% names(trunks)) trunks$sapflow_ref <- NA_real_
if (!"hull_perimeter_cm" %in% names(trunks)) stop("Master layer is missing 'hull_perimeter_cm'.")
if (!"DBH_cm" %in% names(trunks))            stop("Master layer is missing 'DBH_cm'.")
if (!"convexity" %in% names(trunks))         stop("Master layer is missing 'convexity'.")

# ----------------------------
# 3) Compute circular-equivalent sapflow (sapflow_circ) and corrected sapflow_corr
#     - Use r_sapwood if present; otherwise compute local p from sapwood_depth_cm (sapwood mode) or 1 (full_radius)
# ----------------------------
alpha <- u * K
beta  <- (1 - u) * K
mu    <- alpha/(alpha + beta)

have_r_steiner <- "r_steiner" %in% names(trunks)
has_r_sapwood  <- "r_sapwood" %in% names(trunks)
has_swd        <- "sapwood_depth_cm" %in% names(trunks)
have_r_in      <- "r_in" %in% names(trunks)

trunks <- trunks %>%
  rowwise() %>%
  mutate(
    R_cm = as.numeric(DBH_cm) / 2,
    p = if ("r_sapwood" %in% names(cur_data())) {
      as.numeric(r_sapwood)
    } else if (identical(mode, "sapwood")) {
      pmin(pmax(as.numeric(sapwood_depth_cm) / (as.numeric(DBH_cm)/2), 0), 1)
    } else 1,
    
    # --- predictors used for regression ---
    .cp   = pmin(pmax(as.numeric(convexity), eps), 1 - eps),
    .hrel = derive_h_rel(cur_data()),
    
    # === write predicted radii (keep measured ones as-is) ================
    r_in_pred = if (use_regression_for_r) {
      predict_r_in_mu(.cp, ifelse(is.finite(.hrel), .hrel, 0))
    } else as.numeric(r_in),
    r_steiner_pred = if (use_regression_for_r) {
      predict_r_steiner_mu(.cp)
    } else as.numeric(r_steiner),
    
    # physical ordering
    r_steiner_pred = pmin(r_steiner_pred, pmax(0, r_in_pred - 1e-8)),
    
    # --- circular baseline ---
    sapflow_circ = 2 * pi * (R_cm^2) * p * A * (1 - u * p),
    
    # --- Steiner using r_in_pred ---
    sapflow_steiner = 2 * pi * (R_cm^2) * p * A / .cp *
      (1 - .cp * p * mu),
    
    # --- Corrected uses both predicted radii ---
    sapflow_corr = sf_corr_beta(
      R_cm = R_cm, p = p, cp = .cp,
      r_steiner = r_steiner_pred, r_in = r_in_pred,
      alpha = alpha, beta = beta, A = A
    )
  ) %>%
  ungroup() %>%
  select(-any_of(c("p_active", "p", ".cp", ".hrel")))


# ----------------------------
# 4) Bias (unchanged) and write back
# ----------------------------
trunks <- trunks %>%
  mutate(
    SF_circ_bias = ifelse(is.finite(sapflow_ref) & sapflow_ref != 0,
                          (sapflow_circ - sapflow_ref) / sapflow_ref,
                          NA_real_),
    SF_steiner_bias = ifelse(is.finite(sapflow_ref) & sapflow_ref != 0,
                          (sapflow_steiner - sapflow_ref) / sapflow_ref,
                          NA_real_),
    SF_corr_bias = ifelse(is.finite(sapflow_ref) & sapflow_ref != 0,
                          (sapflow_corr - sapflow_ref) / sapflow_ref,
                          NA_real_)
  ) %>%
  # Drop any accidental p_active remnants if they existed in file; do not write local 'p'
  select(-any_of(c("p_active", "p")))

sf::st_write(trunks, master_gpkg, layer = master_layer, delete_layer = TRUE, quiet = F)

message("Updated layer '", master_layer, "' written to: ", master_gpkg)
message("Columns added: sapflow_ref, sapflow_circ, sapflow_steiner, sapflow_corr, ",
        "SF_*_bias, r_in_pred, r_steiner_pred.")



