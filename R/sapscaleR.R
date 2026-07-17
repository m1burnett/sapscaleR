# ----------------------------
# Helpers: robust checks
# ----------------------------
.stopif <- function(cond, msg) if (!isTRUE(cond)) stop(msg, call. = FALSE)
.warnif <- function(cond, msg) if (isTRUE(cond)) warning(msg, call. = FALSE)

clamp01 <- function(x, eps = 1e-9) pmin(pmax(x, eps), 1 - eps)
logit  <- function(x) qlogis(clamp01(x))
invlogit <- function(x) plogis(x)

if (getRversion() >= "2.15.1") {
  utils::globalVariables(c(".x", ".pdf_obs", "tree", "x", "y"))
}

# ----------------------------
# Model coefficients
# ----------------------------
# Point estimates for r_steiner and r_in are from dryad-files/sapflow-estimates.R.
# SEs and precision terms support MC uncertainty draws.
sapscale_models <- list(
  rsteiner = list(
    intercept = -2.885,
    beta_convexity = 0.634,
    se_intercept = 0.205,
    se_beta_convexity = 0.040,
    phi = 20.728,
    phi_se = 7.257
  ),
  rin_convexity_hrel = list(
    intercept = 0.303,
    beta_convexity = 0.297,
    beta_hrel = -1.767,
    se_intercept = 0.293,
    se_beta_convexity = 0.049,
    se_beta_hrel = 0.770,
    phi = 58.474,
    phi_se = 14.376
  ),
  rin_convexity_only = list(
    equation = "Eq. S3.1",
    intercept = -0.306,
    beta_convexity = 0.390,
    se_intercept = 0.131,
    se_beta_convexity = 0.030,
    phi = 51.210,
    phi_se = 12.588
  )
)

# ----------------------------
# Trunk builder
# ----------------------------
# Inputs:
# - Provide convexity directly OR compute from perimeter + hull perimeter
# - Provide diameter or radius. Legacy DBH_cm/Rcirc_cm/Ph_cm names still work.
make_trunk <- function(
    convexity = NULL,
    perimeter = NULL,       # true perimeter
    diameter = NULL,
    radius = NULL,
    hull_perimeter = NULL,  # convex hull perimeter
    perimeter_cm = NULL,
    DBH_cm = NULL,
    Ph_cm = NULL,
    Rcirc_cm = NULL
) {
  if (is.null(perimeter) && !is.null(perimeter_cm)) perimeter <- perimeter_cm
  if (is.null(diameter) && !is.null(DBH_cm)) diameter <- DBH_cm
  if (is.null(radius) && !is.null(Rcirc_cm)) radius <- Rcirc_cm
  if (is.null(hull_perimeter) && !is.null(Ph_cm)) hull_perimeter <- Ph_cm
  
  if (is.null(radius) && !is.null(diameter)) radius <- diameter / 2
  if (is.null(diameter) && !is.null(radius)) diameter <- 2 * radius
  if (is.null(hull_perimeter) && !is.null(diameter)) hull_perimeter <- pi * diameter
  
  # Compute convexity if not provided
  if (is.null(convexity)) {
    .stopif(!is.null(perimeter) && !is.null(hull_perimeter),
            "To compute convexity, provide perimeter and either diameter, radius, or hull_perimeter.")
    convexity <- hull_perimeter / perimeter
  }
  
  .stopif(is.finite(convexity) && convexity > 0 && convexity <= 1,
          "Convexity must be finite and >0 and <=1.")
  .stopif(!is.null(radius) && is.finite(radius) && radius > 0,
          "Either radius or diameter is required, and radius must be positive.")
  
  list(
    convexity = convexity,
    perimeter = perimeter,
    diameter = diameter,
    radius = radius,
    hull_perimeter = hull_perimeter,
    perimeter_cm = perimeter,
    DBH_cm = diameter,
    Ph_cm = hull_perimeter,
    Rcirc_cm = radius
  )
}

make_geom <- make_trunk

# ----------------------------
# Profile objects
# ----------------------------
# A beta curve is a shape on [0, 1]. SF_* eval_mode decides where that
# interval is located in the trunk.
make_beta_curve <- function(mu, K, peak = 1,
                            vcov_log = NULL, se_alpha = NULL, se_beta = NULL) {
  .stopif(is.finite(mu) && mu > 0 && mu < 1, "mu must be in (0,1).")
  .stopif(is.finite(K) && K > 0, "K must be > 0.")
  if (!is.null(vcov_log)) {
    .stopif(is.matrix(vcov_log) && all(dim(vcov_log) == c(2, 2)),
            "vcov_log must be a 2x2 covariance matrix for log(alpha), log(beta).")
  }
  a <- mu * K
  b <- (1 - mu) * K
  list(
    type = "beta_curve",
    mu = mu,
    K = K,
    alpha = a,
    beta = b,
    peak = peak,
    vcov_log = vcov_log,
    se_alpha = se_alpha,
    se_beta = se_beta,
    js_fun = function(x) {
      out <- numeric(length(x))
      ok <- is.finite(x) & x >= 0 & x <= 1
      x01 <- pmin(pmax(x[ok], 1e-9), 1 - 1e-9)
      out[ok] <- peak * dbeta(x01, shape1 = a, shape2 = b)
      out
    }
  )
}

# Legacy wrapper: older code may still create profiles with embedded modes.
make_beta_profile <- function(mode = NULL, mu, K, peak = 1,
                              vcov_log = NULL, se_alpha = NULL, se_beta = NULL,
                              active_max = NULL, active_depth = NULL) {
  .stopif(!is.null(mode), "mode is required: use 'r', 'z', or 'd'.")
  mode <- match.arg(mode, c("r", "z", "d"))
  if (is.null(active_max) && !is.null(active_depth)) active_max <- active_depth
  .stopif(is.finite(mu) && mu > 0 && mu < 1, "mu must be in (0,1).")
  .stopif(is.finite(K) && K > 0, "K must be > 0.")
  if (mode == "d") {
    .stopif(!is.null(active_max) && is.finite(active_max) && active_max > 0,
            "active_depth must be provided and positive for mode='d'.")
  }
  if (!is.null(vcov_log)) {
    .stopif(is.matrix(vcov_log) && all(dim(vcov_log) == c(2, 2)),
            "vcov_log must be a 2x2 covariance matrix for log(alpha), log(beta).")
  }
  a <- mu * K
  b <- (1 - mu) * K
  list(
    mode = mode,
    active_max = if (mode == "d") active_max else 1,
    type = "beta",
    mu = mu,
    K = K,
    alpha = a,
    beta = b,
    peak = peak,
    vcov_log = vcov_log,
    se_alpha = se_alpha,
    se_beta = se_beta,
    js_fun = function(x) {
      x01 <- if (mode == "d") x / active_max else x
      x01 <- pmin(pmax(x01, 0), 1)
      peak * dbeta(x01, shape1 = a, shape2 = b)
    }
  )
}

make_custom_profile <- function(mode = c("r","z","d"), js_fun,
                                active_max = NULL, active_depth = NULL,
                                js_se = NULL, js_se_fun = NULL,
                                draw_min = 0) {
  mode <- match.arg(mode)
  if (is.null(active_max) && !is.null(active_depth)) active_max <- active_depth
  .stopif(is.function(js_fun), "js_fun must be a function(x)->Js.")
  if (!is.null(js_se)) {
    .stopif(length(js_se) == 1 && is.finite(js_se) && js_se >= 0,
            "js_se must be NULL or a non-negative finite scalar.")
  }
  if (!is.null(js_se_fun)) {
    .stopif(is.function(js_se_fun), "js_se_fun must be NULL or a function(x)->SE.")
  }
  if (!is.null(draw_min)) {
    .stopif(length(draw_min) == 1 && is.finite(draw_min),
            "draw_min must be NULL or a finite scalar.")
  }
  list(
    type = "custom",
    mode = mode,
    js_fun = js_fun,
    active_max = active_max,
    js_se = js_se,
    js_se_fun = js_se_fun,
    draw_min = draw_min
  )
}

beta_peak_muK <- function(mu, K) {
  .stopif(is.finite(mu) && mu > 0 && mu < 1, "mu must be in (0,1).")
  .stopif(is.finite(K) && K > 0, "K must be > 0.")
  a <- mu * K
  b <- (1 - mu) * K
  p <- a - 1
  q <- b - 1
  if (!is.finite(p) || !is.finite(q) || p <= 0 || q <= 0) return(NA_real_)
  exp(p * log(p) + q * log(q) - (p + q) * log(p + q) - lbeta(a, b))
}

.beta_unitpeak_pq <- function(x, lp, lq, eps = 1e-6) {
  x <- pmin(pmax(x, eps), 1 - eps)
  p <- exp(lp)
  q <- exp(lq)
  p_log_p <- function(z) ifelse(z <= 1e-12, 0, z * log(z))
  ln_g <- p * log(x) + q * log1p(-x)
  ln_gmax <- p_log_p(p) + p_log_p(q) - (p + q) * log(p + q)
  exp(ln_g - ln_gmax)
}

.fit_unitpeak_beta_optim <- function(x, y, eps = 1e-6, hessian = FALSE) {
  cand <- log(c(0.01, 0.05, 0.1, 0.3, 1, 2, 5))
  grid <- expand.grid(lp = cand, lq = cand)
  rss <- vapply(seq_len(nrow(grid)), function(i) {
    yhat <- .beta_unitpeak_pq(x, grid$lp[i], grid$lq[i], eps = eps)
    if (any(!is.finite(yhat))) return(Inf)
    sum((y - yhat)^2)
  }, numeric(1))
  start <- as.numeric(grid[which.min(rss), c("lp", "lq")])
  opt <- stats::optim(start, function(par) {
    yhat <- .beta_unitpeak_pq(x, par[1], par[2], eps = eps)
    sum((y - yhat)^2)
  }, method = "BFGS", hessian = hessian)
  list(lp = opt$par[1], lq = opt$par[2],
       start_lp = start[1], start_lq = start[2],
       fit = opt)
}

.vcov_log_ab_from_lpq <- function(lp, lq, V_lpq) {
  if (is.null(V_lpq) || !is.matrix(V_lpq) || any(dim(V_lpq) < c(2, 2)) ||
      !all(is.finite(V_lpq))) {
    return(NULL)
  }
  p_tmp <- exp(lp)
  q_tmp <- exp(lq)
  J <- diag(c(p_tmp / (p_tmp + 1), q_tmp / (q_tmp + 1)), 2)
  V <- J %*% V_lpq[seq_len(2), seq_len(2)] %*% J
  if (!all(is.finite(V))) return(NULL)
  V
}

fit_sapflux_beta_profile <- function(
    data,
    scale = c("z", "r", "d"),
    tree_col = "Tree.id",
    sensor_col = "Sensor.id",
    depth_col = "T.depth.cm",
    flux_col = "Avg.SFD.cm.h",
    dbh_col = "DBH.cm",
    dsw_col = NULL,
    r_col = NULL,
    z_col = NULL,
    active_depth = NULL,
    active_depth_cm = NULL,
    normalize = c("sensor_max", "none"),
    exclude_sensor_ids = NULL,
    add_endpoint = TRUE,
    endpoint_x = 0.999,
    endpoint_y = 0.001,
    clamp_y = c(0.001, 0.999),
    method = c("nlme", "optim"),
    eps = 1e-6
) {
  scale <- match.arg(scale)
  normalize <- match.arg(normalize)
  method <- match.arg(method)
  if (is.null(active_depth) && !is.null(active_depth_cm)) active_depth <- active_depth_cm
  
  if (scale == "z" && is.null(z_col) && is.null(dsw_col)) {
    stop("scale='z' requires either z_col or dsw_col.", call. = FALSE)
  }
  if (scale == "r" && is.null(r_col) && is.null(dbh_col)) {
    stop("scale='r' requires either r_col or dbh_col.", call. = FALSE)
  }
  
  need <- c(sensor_col, depth_col, flux_col)
  if (!is.null(tree_col)) need <- c(need, tree_col)
  if (scale == "z" && is.null(z_col)) need <- c(need, dsw_col)
  if (scale == "r" && is.null(r_col)) need <- c(need, dbh_col)
  missing_cols <- setdiff(need, names(data))
  .stopif(length(missing_cols) == 0,
          paste("Missing required column(s):", paste(missing_cols, collapse = ", ")))
  
  df <- data
  if (!is.null(exclude_sensor_ids)) {
    df <- df[!(df[[sensor_col]] %in% exclude_sensor_ids), , drop = FALSE]
  }
  
  sensor <- df[[sensor_col]]
  depth <- as.numeric(df[[depth_col]])
  flux <- as.numeric(df[[flux_col]])
  
  if (normalize == "sensor_max") {
    max_by_sensor <- ave(flux, sensor, FUN = function(v) max(v, na.rm = TRUE))
    y <- flux / max_by_sensor
  } else {
    y <- flux
  }
  y <- pmin(pmax(y, clamp_y[1]), clamp_y[2])
  
  active_depth_values <- NULL
  if (scale == "z") {
    x <- if (!is.null(z_col)) as.numeric(df[[z_col]]) else depth / as.numeric(df[[dsw_col]])
    active_depth_values <- if (!is.null(dsw_col)) as.numeric(df[[dsw_col]]) else NA_real_
  } else if (scale == "r") {
    x <- if (!is.null(r_col)) as.numeric(df[[r_col]]) else depth / (as.numeric(df[[dbh_col]]) / 2)
    active_depth_values <- if (!is.null(dsw_col) && dsw_col %in% names(df)) as.numeric(df[[dsw_col]]) else NA_real_
  } else {
    if (is.null(active_depth)) active_depth <- max(depth, na.rm = TRUE)
    .stopif(is.finite(active_depth) && active_depth > 0,
            "active_depth must be positive for scale='d'.")
    active_depth_values <- rep(active_depth, nrow(df))
    x <- depth / active_depth
  }
  
  fit_df <- data.frame(
    tree = if (!is.null(tree_col)) as.factor(df[[tree_col]]) else factor("all"),
    sensor = sensor,
    x = pmin(pmax(x, eps), 1 - eps),
    y = y,
    depth = depth,
    active_depth = active_depth_values,
    dbh = if (!is.null(dbh_col) && dbh_col %in% names(df)) as.numeric(df[[dbh_col]]) else NA_real_
  )
  fit_df <- fit_df[is.finite(fit_df$x) & is.finite(fit_df$y), , drop = FALSE]
  
  if (isTRUE(add_endpoint)) {
    endpoints <- do.call(rbind, lapply(split(fit_df, fit_df$sensor), function(sdf) {
      row <- sdf[1, , drop = FALSE]
      if (scale == "z" && is.finite(row$active_depth) && is.finite(row$dbh) &&
          row$active_depth / (row$dbh / 2) >= 1) {
        return(NULL)
      }
      row$x <- pmin(pmax(endpoint_x, eps), 1 - eps)
      row$y <- endpoint_y
      if (scale == "r" && is.finite(row$active_depth) && is.finite(row$dbh)) {
        row$x <- pmin(pmax(row$active_depth / (row$dbh / 2), eps), 1 - eps)
      }
      row
    }))
    if (!is.null(endpoints) && nrow(endpoints) > 0) fit_df <- rbind(fit_df, endpoints)
  }
  
  fit_used <- NULL
  if (method == "nlme" && requireNamespace("nlme", quietly = TRUE) &&
      length(unique(fit_df$tree)) > 1) {
    dat <- nlme::groupedData(y ~ x | tree, data = fit_df)
    start_fit <- .fit_unitpeak_beta_optim(fit_df$x, fit_df$y, eps = eps)
    start <- c(lp = start_fit$start_lp, lq = start_fit$start_lq)
    nlme_formula <- y ~ .beta_unitpeak_pq(x, lp, lq)
    environment(nlme_formula) <- environment(.beta_unitpeak_pq)
    fit_used <- try(nlme::nlme(
      nlme_formula,
      data = dat,
      fixed = lp + lq ~ 1,
      random = nlme::pdDiag(lq ~ 1),
      groups = ~ tree,
      start = start,
      control = nlme::nlmeControl(
        pnlsTol = 0.1,
        pnlsMaxIter = 50,
        msMaxIter = 300,
        returnObject = TRUE
      )
    ), silent = TRUE)
  }
  
  if (!inherits(fit_used, "nlme")) {
    opt_fit <- .fit_unitpeak_beta_optim(fit_df$x, fit_df$y, eps = eps, hessian = TRUE)
    lp <- opt_fit$lp
    lq <- opt_fit$lq
    method_used <- "optim"
    fit_used <- opt_fit$fit
    vcov_log <- NULL
    se_alpha <- NULL
    se_beta <- NULL
    H <- fit_used$hessian
    if (is.matrix(H) && all(dim(H) == c(2, 2)) && all(is.finite(H))) {
      V_lpq <- try(2 * (fit_used$value / max(nrow(fit_df) - 2, 1)) * solve(H),
                   silent = TRUE)
      if (!inherits(V_lpq, "try-error")) {
        vcov_log <- .vcov_log_ab_from_lpq(lp, lq, V_lpq)
        if (!is.null(vcov_log)) {
          se_log <- sqrt(pmax(diag(vcov_log), 0))
          se_alpha <- (exp(lp) + 1) * se_log[1]
          se_beta <- (exp(lq) + 1) * se_log[2]
        }
      }
    }
  } else {
    fx <- nlme::fixef(fit_used)
    lp <- unname(fx["lp"])
    lq <- unname(fx["lq"])
    method_used <- "nlme"
    vcov_log <- NULL
    se_alpha <- NULL
    se_beta <- NULL
    V_lpq <- try({
      V <- if (!is.null(fit_used$varFix)) {
        as.matrix(fit_used$varFix)
      } else {
        as.matrix(stats::vcov(fit_used))
      }
      if (all(c("lp", "lq") %in% rownames(V)) &&
          all(c("lp", "lq") %in% colnames(V))) {
        V[c("lp", "lq"), c("lp", "lq")]
      } else {
        V[seq_len(2), seq_len(2)]
      }
    }, silent = TRUE)
    if (!inherits(V_lpq, "try-error") && all(is.finite(V_lpq))) {
      vcov_log <- .vcov_log_ab_from_lpq(lp, lq, V_lpq)
      if (!is.null(vcov_log)) {
        se_log <- sqrt(pmax(diag(vcov_log), 0))
        se_alpha <- (exp(lp) + 1) * se_log[1]
        se_beta <- (exp(lq) + 1) * se_log[2]
      }
    }
  }
  
  p <- exp(lp)
  q <- exp(lq)
  K <- p + q + 2
  mu <- (p + 1) / K
  fmax <- beta_peak_muK(mu, K)
  active_max <- if (scale == "d") active_depth else NULL
  
  list(
    scale = scale,
    mu = mu,
    K = K,
    alpha = mu * K,
    beta = (1 - mu) * K,
    fmax = fmax,
    unit_peak_multiplier = if (is.finite(fmax)) 1 / fmax else NA_real_,
    vcov_log = vcov_log,
    se_alpha = se_alpha,
    se_beta = se_beta,
    method = method_used,
    fit = fit_used,
    fit_data = fit_df,
    profile = make_beta_curve(
      mu = mu, K = K, peak = 1,
      vcov_log = vcov_log, se_alpha = se_alpha, se_beta = se_beta
    ),
    unit_peak_profile = make_beta_curve(
      mu = mu,
      K = K,
      peak = if (is.finite(fmax)) 1 / fmax else 1,
      vcov_log = vcov_log,
      se_alpha = se_alpha,
      se_beta = se_beta
    )
  )
}

# ----------------------------
# Depth mapping
# ----------------------------
# Convert evaluation grid in chosen "eval_mode" into r and into the profile's native mode.
# Inputs:
# - geom: from make_trunk()
# - dsw: required when profile mode is 'z' (or if you ever evaluate on z),
#        since p=dsw/Rcirc and r = p*z.
.resolve_dsw <- function(dsw = NULL, dsw_cm = NULL) {
  if (is.null(dsw) && !is.null(dsw_cm)) dsw <- dsw_cm
  dsw
}

.resolve_dsw_se <- function(dsw_se = NULL, dsw_se_cm = NULL) {
  if (is.null(dsw_se) && !is.null(dsw_se_cm)) dsw_se <- dsw_se_cm
  dsw_se
}

map_depths <- function(x, eval_mode = c("r","z"), profile, geom,
                       dsw = NULL, dsw_cm = NULL) {
  eval_mode <- match.arg(eval_mode)
  prof_mode <- profile$mode
  dsw <- .resolve_dsw(dsw = dsw, dsw_cm = dsw_cm)
  
  Rcirc <- geom$Rcirc_cm
  .stopif(is.finite(Rcirc) && Rcirc > 0, "trunk must include radius.")
  
  need_dsw <- (eval_mode == "z") || (prof_mode == "z")
  if (need_dsw) .stopif(!is.null(dsw) && is.finite(dsw) && dsw > 0,
                        "dsw is required when using z (either eval_mode='z' or profile mode='z').")
  
  p <- if (need_dsw) pmin(dsw / Rcirc, 1) else NA_real_
  
  # Convert eval grid to r
  r <- switch(
    eval_mode,
    r = x,
    z = p * x
  )
  
  # Convert eval grid to profile-native coordinate
  x_prof <- switch(
    prof_mode,
    r = r,
    z = if (eval_mode == "z") x else r / p,
    d = r * Rcirc
  )
  
  list(r = r, x_prof = x_prof, p = p)
}

map_profile_depths <- function(r, eval_mode = c("r", "z", "d"), profile, geom,
                               dsw = NULL, dsw_cm = NULL,
                               active_depth = NULL) {
  eval_mode <- match.arg(eval_mode)
  dsw <- .resolve_dsw(dsw = dsw, dsw_cm = dsw_cm)
  Rcirc <- geom$Rcirc_cm
  .stopif(is.finite(Rcirc) && Rcirc > 0, "trunk must include radius.")
  
  if (identical(profile$type, "beta_curve")) {
    p <- NA_real_
    x_prof <- switch(
      eval_mode,
      r = r,
      z = {
        .stopif(!is.null(dsw) && is.finite(dsw) && dsw > 0,
                "dsw is required when eval_mode='z'.")
        p <- pmin(dsw / Rcirc, 1)
        r / p
      },
      d = {
        .stopif(!is.null(active_depth) && is.finite(active_depth) && active_depth > 0,
                "active_depth is required when eval_mode='d'.")
        r * Rcirc / active_depth
      }
    )
    active_r_max <- switch(
      eval_mode,
      r = 1,
      z = p,
      d = min(active_depth / Rcirc, 1)
    )
    return(list(r = r, x_prof = x_prof, p = p, active_r_max = active_r_max))
  }
  
  # Legacy mode-aware profiles, including make_custom_profile().
  mp <- map_depths(r, eval_mode = "r", profile = profile, geom = geom, dsw = dsw)
  active_r_max <- if (is.finite(mp$p)) mp$p else NULL
  if (is.null(active_r_max) && identical(profile$mode, "d") &&
      !is.null(profile$active_max) && is.finite(profile$active_max)) {
    active_r_max <- min(profile$active_max / Rcirc, 1)
  }
  mp$active_r_max <- active_r_max
  mp
}

# ----------------------------
# rsteiner / rin prediction hooks
# ----------------------------
# Provide your fitted model coefficients here.
# rsteiner model in draft is a beta regression vs logit(convexity).
predict_rsteiner <- function(convexity, models = sapscale_models) {
  .stopif(!is.null(models$rsteiner), "models$rsteiner must be provided.")
  m <- models$rsteiner
  eta <- m$intercept + m$beta_convexity * logit(convexity)
  invlogit(eta)
}

# rin can use convexity-only or convexity+hrel; warn if corr is requested but h_rel is missing.
.has_hrel <- function(h_rel) {
  !is.null(h_rel) && length(h_rel) == 1 && is.finite(h_rel)
}

.select_rin_model <- function(h_rel = NULL, models = sapscale_models,
                              prefer_hrel = TRUE, warn_missing_hrel = TRUE) {
  if (prefer_hrel && .has_hrel(h_rel) && !is.null(models$rin_convexity_hrel)) {
    return(list(model = models$rin_convexity_hrel, uses_hrel = TRUE))
  }

  .warnif(prefer_hrel && !.has_hrel(h_rel) && warn_missing_hrel,
          "h_rel not provided; using convexity-only r_in model (Eq. S3.1).")
  .stopif(!is.null(models$rin_convexity_only),
          "models$rin_convexity_only must be provided for convexity-only rin fallback.")
  list(model = models$rin_convexity_only, uses_hrel = FALSE)
}

predict_rin <- function(convexity, h_rel = NULL, models = sapscale_models,
                        prefer_hrel = TRUE, warn_missing_hrel = TRUE) {
  choice <- .select_rin_model(
    h_rel = h_rel,
    models = models,
    prefer_hrel = prefer_hrel,
    warn_missing_hrel = warn_missing_hrel
  )
  m <- choice$model
  eta <- m$intercept + m$beta_convexity * logit(convexity)
  if (choice$uses_hrel) eta <- eta + m$beta_hrel * h_rel
  invlogit(eta)
}

# Convenience: compute h_rel from h_max and diameter
hrel_from_hmax <- function(h_max = NULL, diameter = NULL,
                           h_max_cm = NULL, DBH_cm = NULL) {
  if (is.null(h_max) && !is.null(h_max_cm)) h_max <- h_max_cm
  if (is.null(diameter) && !is.null(DBH_cm)) diameter <- DBH_cm
  .stopif(is.finite(h_max) && h_max >= 0, "h_max must be >= 0.")
  .stopif(is.finite(diameter) && diameter > 0, "diameter must be > 0.")
  h_max / diameter
}

# ----------------------------
# Core integrators for SFcirc / SFsteiner / SFcorr
# ----------------------------
# Numerical integration over r in [0,1] with a supplied "g(r)" using trapezoids.
integrate_sf <- function(r_grid, Js_grid, g_grid, scale_area = 1) {
  o <- order(r_grid)
  r <- r_grid[o]; Js <- Js_grid[o]; g <- g_grid[o]
  dr <- diff(r)
  mid <- (Js[-1] * g[-1] + Js[-length(Js)] * g[-length(g)]) / 2
  scale_area * sum(mid * dr)
}

# g for circle in r-space
g_circ <- function(r) pmax(1 - r, 0)

# Shallow-erosion / Steiner-valid g(r) hook
# NOTE: Keep this as the simple placeholder until you drop in your finalized expression.
g_steiner <- function(r, convexity) pmax((1 / convexity) - r, 0)

# Piecewise linear g_corr(r): from g(rsteiner) to 0 at rin.
# active_r_max is accepted for backwards compatibility but does not move the
# geometric taper endpoint; sap flux profiles define the active sapwood domain.
g_corr <- function(r, convexity, rsteiner, rin, active_r_max = NULL) {
  .stopif(rin > 0 && rin <= 1, "rin must be in (0,1].")
  .stopif(rsteiner > 0 && rsteiner < rin, "Need 0<rsteiner<rin for SFcorr.")
  
  g_rs <- g_steiner(rsteiner, convexity)
  
  out <- numeric(length(r))
  for (i in seq_along(r)) {
    if (r[i] <= rsteiner) {
      out[i] <- g_steiner(r[i], convexity)
    } else if (r[i] <= rin) {
      out[i] <- g_rs * (1 - (r[i] - rsteiner) / (rin - rsteiner))
    } else {
      out[i] <- 0
    }
  }
  pmax(out, 0)
}

# ----------------------------
# Main public functions: SFcirc, SFsteiner, SFcorr
# ----------------------------
.resolve_h_rel <- function(h_rel = NULL, h_max = NULL, h_max_cm = NULL, geom) {
  if (is.null(h_max) && !is.null(h_max_cm)) h_max <- h_max_cm
  if (is.null(h_rel) && !is.null(h_max)) {
    .stopif(!is.null(geom$diameter) && is.finite(geom$diameter) && geom$diameter > 0,
            "To compute h_rel from h_max, geom must include diameter.")
    h_rel <- hrel_from_hmax(h_max = h_max, diameter = geom$diameter)
  }
  h_rel
}

.SF_circ_value <- function(profile, geom, eval_mode = c("r","z","d"),
                           dsw = NULL, dsw_cm = NULL,
                           active_depth = NULL, n = 2001) {
  eval_mode <- match.arg(eval_mode)
  dsw <- .resolve_dsw(dsw = dsw, dsw_cm = dsw_cm)
  Rcirc <- geom$Rcirc_cm
  
  # Evaluate on r-grid; profile can be in r or z via map_depths()
  r_grid <- seq(0, 1, length.out = n)
  mp <- map_profile_depths(r_grid, eval_mode = eval_mode, profile = profile,
                           geom = geom, dsw = dsw, active_depth = active_depth)
  Js <- profile$js_fun(mp$x_prof)
  
  # Area scale for the r-integral form
  scale_area <- 2 * pi * Rcirc^2
  
  integrate_sf(r_grid, Js, g_circ(r_grid), scale_area = scale_area)
}

.SF_steiner_value <- function(profile, geom, eval_mode = c("r","z","d"),
                              dsw = NULL, dsw_cm = NULL,
                              active_depth = NULL, n = 2001) {
  eval_mode <- match.arg(eval_mode)
  dsw <- .resolve_dsw(dsw = dsw, dsw_cm = dsw_cm)
  Rcirc <- geom$Rcirc_cm
  convexity <- geom$convexity
  
  r_grid <- seq(0, 1, length.out = n)
  mp <- map_profile_depths(r_grid, eval_mode = eval_mode, profile = profile,
                           geom = geom, dsw = dsw, active_depth = active_depth)
  Js <- profile$js_fun(mp$x_prof)
  
  scale_area <- 2 * pi * Rcirc^2
  
  integrate_sf(r_grid, Js, g_steiner(r_grid, convexity), scale_area = scale_area)
}

.SF_corr_value <- function(profile, geom, models = sapscale_models,
                           eval_mode = c("r","z","d"),
                           dsw = NULL, dsw_cm = NULL,
                           active_depth = NULL,
                           h_rel = NULL, h_max = NULL, h_max_cm = NULL,
                           prefer_hrel = TRUE,
                           n = 2001,
                           rsteiner = NULL, rin = NULL,
                           warn_missing_hrel = TRUE) {
  eval_mode <- match.arg(eval_mode)
  dsw <- .resolve_dsw(dsw = dsw, dsw_cm = dsw_cm)
  convexity <- geom$convexity
  Rcirc <- geom$Rcirc_cm
  
  h_rel <- .resolve_h_rel(h_rel = h_rel, h_max = h_max, h_max_cm = h_max_cm, geom = geom)
  
  if (is.null(rsteiner)) rsteiner <- predict_rsteiner(convexity, models)
  if (is.null(rin)) {
    rin <- predict_rin(
      convexity,
      h_rel = h_rel,
      models = models,
      prefer_hrel = prefer_hrel,
      warn_missing_hrel = warn_missing_hrel
    )
  }
  
  r_grid <- seq(0, 1, length.out = n)
  mp <- map_profile_depths(r_grid, eval_mode = eval_mode, profile = profile,
                           geom = geom, dsw = dsw, active_depth = active_depth)
  Js <- profile$js_fun(mp$x_prof)
  active_r_max <- mp$active_r_max
  
  scale_area <- 2 * pi * Rcirc^2
  
  integrate_sf(r_grid, Js, g_corr(r_grid, convexity, rsteiner, rin,
                                  active_r_max = active_r_max),
               scale_area = scale_area)
}

.check_mc_args <- function(n_draws, conf) {
  .stopif(length(n_draws) == 1 && is.finite(n_draws) && n_draws > 0,
          "n_draws must be a positive number.")
  .stopif(length(conf) == 1 && is.finite(conf) && conf > 0 && conf < 1,
          "conf must be a number between 0 and 1.")
  as.integer(n_draws)
}

.mc_summary <- function(draws, n_draws, conf, keep_draws,
                        uncertainty_sources = character(0)) {
  valid <- is.finite(draws)
  good <- draws[valid]
  if (!length(good)) stop("All MC draws were invalid; check inputs.", call. = FALSE)
  
  alpha <- (1 - conf) / 2
  ci <- stats::quantile(good, probs = c(alpha, 1 - alpha), names = FALSE)
  names(ci) <- c("lower", "upper")
  
  out <- list(
    estimate = mean(good),
    se = stats::sd(good),
    ci = ci,
    conf = conf,
    uncertainty_sources = uncertainty_sources,
    n_draws = n_draws,
    n_valid = length(good),
    n_invalid = n_draws - length(good)
  )
  if (isTRUE(keep_draws)) out$draws <- good
  out
}

.is_pos_number <- function(x) {
  !is.null(x) && length(x) == 1 && is.finite(x) && x > 0
}

.profile_value_vector <- function(value, x, name) {
  value <- as.numeric(value)
  if (length(value) == 1 && length(x) != 1) value <- rep(value, length(x))
  .stopif(length(value) == length(x),
          paste0(name, " must return a scalar or a vector the same length as x."))
  value
}

.profile_vcov_log <- function(profile) {
  if (!is.null(profile$vcov_log) &&
      is.matrix(profile$vcov_log) &&
      all(dim(profile$vcov_log) == c(2, 2)) &&
      all(is.finite(profile$vcov_log))) {
    return(profile$vcov_log)
  }
  nms <- names(profile)
  if (all(c("vcov_log_11", "vcov_log_12", "vcov_log_22") %in% nms) &&
      all(is.finite(c(profile$vcov_log_11, profile$vcov_log_12, profile$vcov_log_22)))) {
    return(matrix(c(profile$vcov_log_11, profile$vcov_log_12,
                    profile$vcov_log_12, profile$vcov_log_22), 2, 2))
  }
  NULL
}

.has_beta_profile_uncertainty <- function(profile) {
  !is.null(.profile_vcov_log(profile)) ||
    (.is_pos_number(profile$se_alpha) && .is_pos_number(profile$se_beta))
}

.has_custom_profile_uncertainty <- function(profile) {
  is.function(profile$js_se_fun) || .is_pos_number(profile$js_se)
}

.has_profile_uncertainty <- function(profile) {
  .has_beta_profile_uncertainty(profile) || .has_custom_profile_uncertainty(profile)
}

.custom_profile_se <- function(profile, x) {
  if (is.function(profile$js_se_fun)) {
    se <- profile$js_se_fun(x)
  } else if (!is.null(profile$js_se)) {
    se <- profile$js_se
  } else {
    se <- 0
  }
  se <- .profile_value_vector(se, x, "js_se_fun")
  pmax(se, 0)
}

.draw_profile_params <- function(profile, n_draws, include_profile_uncertainty = TRUE) {
  if (!all(c("alpha", "beta") %in% names(profile)) ||
      !.is_pos_number(profile$alpha) || !.is_pos_number(profile$beta)) {
    return(list(alpha = NULL, beta = NULL))
  }
  
  alpha_hat <- profile$alpha
  beta_hat <- profile$beta
  
  if (isTRUE(include_profile_uncertainty)) {
    V <- .profile_vcov_log(profile)
    if (!is.null(V)) {
      mu <- c(log(alpha_hat), log(beta_hat))
      L <- try(chol(V), silent = TRUE)
      if (inherits(L, "try-error")) {
        V <- diag(pmax(diag(V), 0), 2)
        L <- chol(V + diag(1e-12, 2))
      }
      Z <- matrix(stats::rnorm(n_draws * 2), ncol = 2)
      theta <- sweep(Z %*% t(L), 2, mu, "+")
      return(list(alpha = exp(theta[, 1]), beta = exp(theta[, 2])))
    }
    
    if (.is_pos_number(profile$se_alpha) && .is_pos_number(profile$se_beta)) {
      sd_loga <- profile$se_alpha / alpha_hat
      sd_logb <- profile$se_beta / beta_hat
      return(list(
        alpha = exp(stats::rnorm(n_draws, mean = log(alpha_hat), sd = sd_loga)),
        beta = exp(stats::rnorm(n_draws, mean = log(beta_hat), sd = sd_logb))
      ))
    }
  }
  
  list(alpha = rep(alpha_hat, n_draws), beta = rep(beta_hat, n_draws))
}

.make_beta_draw_profile <- function(profile, alpha, beta) {
  peak <- if (is.null(profile$peak)) 1 else profile$peak
  if (identical(profile$type, "beta_curve")) {
    return(make_beta_curve(mu = alpha / (alpha + beta),
                           K = alpha + beta,
                           peak = peak))
  }
  active_max <- if (is.null(profile$active_max)) NULL else profile$active_max
  make_custom_profile(
    mode = profile$mode,
    active_max = active_max,
    js_fun = function(x) {
      x01 <- if (identical(profile$mode, "d")) x / active_max else x
      x01 <- pmin(pmax(x01, 0), 1)
      peak * stats::dbeta(x01, shape1 = alpha, shape2 = beta)
    }
  )
}

.make_custom_draw_profile <- function(profile) {
  out <- profile
  base_fun <- profile$js_fun
  draw_min <- profile$draw_min
  out$js_fun <- function(x) {
    mu <- .profile_value_vector(base_fun(x), x, "js_fun")
    se <- .custom_profile_se(profile, x)
    draw <- stats::rnorm(length(mu), mean = mu, sd = se)
    if (!is.null(draw_min)) draw <- pmax(draw, draw_min)
    draw
  }
  out
}

.draw_dsw <- function(dsw_cm, dsw_se_cm, n_draws, eval_mode) {
  if (eval_mode != "z" || !.is_pos_number(dsw_se_cm)) return(dsw_cm)
  .stopif(.is_pos_number(dsw_cm),
          "dsw must be provided and positive when drawing dsw uncertainty.")
  
  dsw_draw <- stats::rnorm(n_draws, mean = dsw_cm, sd = dsw_se_cm)
  bad <- which(!is.finite(dsw_draw) | dsw_draw <= 0)
  iter <- 0
  while (length(bad) > 0 && iter < 50) {
    dsw_draw[bad] <- stats::rnorm(length(bad), mean = dsw_cm, sd = dsw_se_cm)
    bad <- which(!is.finite(dsw_draw) | dsw_draw <= 0)
    iter <- iter + 1
  }
  dsw_draw[!is.finite(dsw_draw) | dsw_draw <= 0] <- min(dsw_cm, 1e-6)
  dsw_draw
}

.draw_normal_param <- function(point, se, n_draws) {
  .stopif(!is.null(point) && length(point) == 1 && is.finite(point),
          "Model coefficient point estimates must be finite.")
  if (!.is_pos_number(se)) return(rep(point, n_draws))
  stats::rnorm(n_draws, mean = point, sd = se)
}

.draw_positive_param <- function(point, se, n_draws) {
  if (!.is_pos_number(point)) return(rep(NA_real_, n_draws))
  if (!.is_pos_number(se)) return(rep(point, n_draws))
  sdlog <- sqrt(log1p((se / point)^2))
  meanlog <- log(point) - 0.5 * sdlog^2
  stats::rlnorm(n_draws, meanlog = meanlog, sdlog = sdlog)
}

.draw_rmodel <- function(model, convexity, h_rel = NULL, n_draws,
                         include_residual = TRUE) {
  b0 <- .draw_normal_param(model$intercept, model$se_intercept, n_draws)
  b1 <- .draw_normal_param(model$beta_convexity, model$se_beta_convexity, n_draws)
  eta <- b0 + b1 * logit(convexity)
  
  if (!is.null(model$beta_hrel) && .has_hrel(h_rel)) {
    b2 <- .draw_normal_param(model$beta_hrel, model$se_beta_hrel, n_draws)
    eta <- eta + b2 * h_rel
  }
  
  mu <- clamp01(invlogit(eta))
  phi <- .draw_positive_param(model$phi, model$phi_se, n_draws)
  
  if (isTRUE(include_residual) && all(is.finite(phi))) {
    return(stats::rbeta(n_draws, shape1 = mu * phi, shape2 = (1 - mu) * phi))
  }
  mu
}

.check_radius_se <- function(se, name) {
  if (is.null(se)) return(NULL)
  .stopif(length(se) == 1 && is.finite(se) && se >= 0,
          paste0(name, " must be NULL or a non-negative finite scalar."))
  se
}

.draw_user_radius <- function(point, se, n_draws) {
  if (.is_pos_number(se)) {
    return(stats::rnorm(n_draws, mean = point, sd = se))
  }
  rep(point, n_draws)
}

.draw_corr_radii <- function(convexity, h_rel, models, prefer_hrel,
                             n_draws, include_rmodel_residual,
                             rsteiner = NULL, rin = NULL,
                             rsteiner_se = NULL, rin_se = NULL,
                             warn_missing_hrel = TRUE, max_tries = 20) {
  have_rsteiner <- !is.null(rsteiner)
  have_rin <- !is.null(rin)
  rsteiner_se <- .check_radius_se(rsteiner_se, "r_steiner_se")
  rin_se <- .check_radius_se(rin_se, "r_in_se")
  .stopif(is.null(rsteiner_se) || have_rsteiner,
          "r_steiner_se requires a user-supplied r_steiner value.")
  .stopif(is.null(rin_se) || have_rin,
          "r_in_se requires a user-supplied r_in value.")
  if (have_rsteiner) .stopif(length(rsteiner) == 1 && is.finite(rsteiner),
                             "r_steiner must be a finite scalar.")
  if (have_rin) .stopif(length(rin) == 1 && is.finite(rin),
                        "r_in must be a finite scalar.")
  if (have_rsteiner && have_rin &&
      !.is_pos_number(rsteiner_se) && !.is_pos_number(rin_se)) {
    .stopif(rsteiner > 0 && rsteiner < rin && rin <= 1,
            "Need 0 < r_steiner < r_in <= 1.")
    return(list(
      rsteiner = rep(rsteiner, n_draws),
      rin = rep(rin, n_draws)
    ))
  }
  
  if (!have_rsteiner) .stopif(!is.null(models$rsteiner), "models$rsteiner must be provided.")
  rin_choice <- NULL
  if (!have_rin) {
    rin_choice <- .select_rin_model(
      h_rel = h_rel,
      models = models,
      prefer_hrel = prefer_hrel,
      warn_missing_hrel = warn_missing_hrel
    )
  }
  
  rsteiner_draw <- rep(NA_real_, n_draws)
  rin_draw <- rep(NA_real_, n_draws)
  pending <- seq_len(n_draws)
  iter <- 0
  
  while (length(pending) > 0 && iter < max_tries) {
    n_pending <- length(pending)
    rs <- if (have_rsteiner) {
      .draw_user_radius(rsteiner, rsteiner_se, n_pending)
    } else {
      .draw_rmodel(
        models$rsteiner,
        convexity = convexity,
        n_draws = n_pending,
        include_residual = include_rmodel_residual
      )
    }
    ri <- if (have_rin) {
      .draw_user_radius(rin, rin_se, n_pending)
    } else {
      .draw_rmodel(
        rin_choice$model,
        convexity = convexity,
        h_rel = if (rin_choice$uses_hrel) h_rel else NULL,
        n_draws = n_pending,
        include_residual = include_rmodel_residual
      )
    }
    
    ok <- is.finite(rs) & is.finite(ri) & rs > 0 & rs < ri & ri <= 1
    if (any(ok)) {
      idx <- pending[ok]
      rsteiner_draw[idx] <- rs[ok]
      rin_draw[idx] <- ri[ok]
    }
    pending <- pending[!ok]
    iter <- iter + 1
  }
  
  list(rsteiner = rsteiner_draw, rin = rin_draw)
}

.SF_mc <- function(method, profile, geom, models = sapscale_models,
                   eval_mode = c("r","z","d"), dsw = NULL, dsw_cm = NULL,
                   active_depth = NULL,
                   h_rel = NULL, h_max = NULL, h_max_cm = NULL, prefer_hrel = TRUE,
                   rsteiner = NULL, rin = NULL,
                   rsteiner_se = NULL, rin_se = NULL,
                   n = 2001, n_draws = 2000, conf = 0.95,
                   seed = NULL, keep_draws = FALSE, dsw_se = NULL, dsw_se_cm = NULL,
                   include_profile_uncertainty = TRUE,
                   include_rmodel_residual = TRUE) {
  method <- match.arg(method, c("circ", "steiner", "corr"))
  eval_mode <- match.arg(eval_mode)
  dsw <- .resolve_dsw(dsw = dsw, dsw_cm = dsw_cm)
  dsw_se <- .resolve_dsw_se(dsw_se = dsw_se, dsw_se_cm = dsw_se_cm)
  n_draws <- .check_mc_args(n_draws, conf)
  if (!is.null(seed)) set.seed(seed)
  
  profile_draws <- .draw_profile_params(
    profile,
    n_draws = n_draws,
    include_profile_uncertainty = include_profile_uncertainty
  )
  dsw_draw <- .draw_dsw(
    dsw_cm = dsw,
    dsw_se_cm = dsw_se,
    n_draws = n_draws,
    eval_mode = eval_mode
  )
  
  h_rel <- .resolve_h_rel(h_rel = h_rel, h_max = h_max, h_max_cm = h_max_cm, geom = geom)
  profile_uncertainty_active <- isTRUE(include_profile_uncertainty) &&
    .has_profile_uncertainty(profile)
  uncertainty_sources <- character(0)
  if (profile_uncertainty_active) uncertainty_sources <- c(uncertainty_sources, "profile")
  if (eval_mode == "z" && .is_pos_number(dsw_se)) {
    uncertainty_sources <- c(uncertainty_sources, "dsw")
  }
  if (method == "corr") {
    if (is.null(rsteiner) || .is_pos_number(rsteiner_se)) {
      uncertainty_sources <- c(uncertainty_sources, "r_steiner")
    }
    if (is.null(rin) || .is_pos_number(rin_se)) {
      uncertainty_sources <- c(uncertainty_sources, "r_in")
    }
  }
  uncertainty_sources <- unique(uncertainty_sources)

  radii <- NULL
  if (method == "corr") {
    radii <- .draw_corr_radii(
      convexity = geom$convexity,
      h_rel = h_rel,
      models = models,
      prefer_hrel = prefer_hrel,
      n_draws = n_draws,
      include_rmodel_residual = include_rmodel_residual,
      rsteiner = rsteiner,
      rin = rin,
      rsteiner_se = rsteiner_se,
      rin_se = rin_se,
      warn_missing_hrel = TRUE
    )
  }
  
  draws <- rep(NA_real_, n_draws)
  for (i in seq_len(n_draws)) {
    prof_i <- profile
    if (!is.null(profile_draws$alpha)) {
      prof_i <- .make_beta_draw_profile(profile, profile_draws$alpha[i], profile_draws$beta[i])
    } else if (profile_uncertainty_active && .has_custom_profile_uncertainty(profile)) {
      prof_i <- .make_custom_draw_profile(profile)
    }
    dsw_i <- if (is.null(dsw_draw)) NULL else dsw_draw[min(i, length(dsw_draw))]
    
    val <- try({
      if (method == "circ") {
        .SF_circ_value(profile = prof_i, geom = geom, eval_mode = eval_mode,
                       dsw = dsw_i, active_depth = active_depth, n = n)
      } else if (method == "steiner") {
        .SF_steiner_value(profile = prof_i, geom = geom, eval_mode = eval_mode,
                          dsw = dsw_i, active_depth = active_depth, n = n)
      } else {
        if (!is.finite(radii$rsteiner[i]) || !is.finite(radii$rin[i])) {
          NA_real_
        } else {
          .SF_corr_value(profile = prof_i, geom = geom, models = models,
                         eval_mode = eval_mode, dsw = dsw_i,
                         active_depth = active_depth,
                         h_rel = h_rel, h_max = NULL, h_max_cm = NULL,
                         prefer_hrel = prefer_hrel, n = n,
                         rsteiner = radii$rsteiner[i], rin = radii$rin[i],
                         warn_missing_hrel = FALSE)
        }
      }
    }, silent = TRUE)
    
    if (!inherits(val, "try-error") && length(val) == 1 && is.finite(val)) {
      draws[i] <- val
    }
  }
  
  .mc_summary(
    draws,
    n_draws = n_draws,
    conf = conf,
    keep_draws = keep_draws,
    uncertainty_sources = uncertainty_sources
  )
}

SF_circ <- function(profile, geom, eval_mode = c("r","z","d"),
                    dsw = NULL, dsw_cm = NULL,
                    active_depth = NULL,
                    n = 2001,
                    uncertainty = c("none", "mc"),
                    n_draws = 2000, conf = 0.95, seed = NULL,
                    keep_draws = FALSE, dsw_se = NULL, dsw_se_cm = NULL,
                    include_profile_uncertainty = TRUE) {
  uncertainty <- match.arg(uncertainty)
  dsw <- .resolve_dsw(dsw = dsw, dsw_cm = dsw_cm)
  dsw_se <- .resolve_dsw_se(dsw_se = dsw_se, dsw_se_cm = dsw_se_cm)
  if (uncertainty == "mc") {
    return(.SF_mc(
      method = "circ",
      profile = profile,
      geom = geom,
      eval_mode = eval_mode,
      dsw = dsw,
      active_depth = active_depth,
      n = n,
      n_draws = n_draws,
      conf = conf,
      seed = seed,
      keep_draws = keep_draws,
      dsw_se = dsw_se,
      include_profile_uncertainty = include_profile_uncertainty
    ))
  }
  .SF_circ_value(profile = profile, geom = geom, eval_mode = eval_mode,
                 dsw = dsw, active_depth = active_depth, n = n)
}

SF_steiner <- function(profile, geom, eval_mode = c("r","z","d"),
                       dsw = NULL, dsw_cm = NULL,
                       active_depth = NULL,
                       n = 2001,
                       uncertainty = c("none", "mc"),
                       n_draws = 2000, conf = 0.95, seed = NULL,
                       keep_draws = FALSE, dsw_se = NULL, dsw_se_cm = NULL,
                       include_profile_uncertainty = TRUE) {
  uncertainty <- match.arg(uncertainty)
  dsw <- .resolve_dsw(dsw = dsw, dsw_cm = dsw_cm)
  dsw_se <- .resolve_dsw_se(dsw_se = dsw_se, dsw_se_cm = dsw_se_cm)
  if (uncertainty == "mc") {
    return(.SF_mc(
      method = "steiner",
      profile = profile,
      geom = geom,
      eval_mode = eval_mode,
      dsw = dsw,
      active_depth = active_depth,
      n = n,
      n_draws = n_draws,
      conf = conf,
      seed = seed,
      keep_draws = keep_draws,
      dsw_se = dsw_se,
      include_profile_uncertainty = include_profile_uncertainty
    ))
  }
  .SF_steiner_value(profile = profile, geom = geom, eval_mode = eval_mode,
                    dsw = dsw, active_depth = active_depth, n = n)
}

SF_corr <- function(profile, geom, models = sapscale_models,
                    eval_mode = c("r","z","d"),
                    dsw = NULL, dsw_cm = NULL,
                    active_depth = NULL,
                    h_rel = NULL, h_max = NULL, h_max_cm = NULL,
                    r_steiner = NULL, r_in = NULL,
                    r_steiner_se = NULL, r_in_se = NULL,
                    prefer_hrel = TRUE,
                    n = 2001,
                    uncertainty = c("none", "mc"),
                    n_draws = 2000, conf = 0.95, seed = NULL,
                    keep_draws = FALSE, dsw_se = NULL, dsw_se_cm = NULL,
                    include_profile_uncertainty = TRUE,
                    include_rmodel_residual = TRUE) {
  uncertainty <- match.arg(uncertainty)
  dsw <- .resolve_dsw(dsw = dsw, dsw_cm = dsw_cm)
  dsw_se <- .resolve_dsw_se(dsw_se = dsw_se, dsw_se_cm = dsw_se_cm)
  if (uncertainty == "mc") {
    return(.SF_mc(
      method = "corr",
      profile = profile,
      geom = geom,
      models = models,
      eval_mode = eval_mode,
      dsw = dsw,
      active_depth = active_depth,
      h_rel = h_rel,
      h_max = h_max,
      h_max_cm = h_max_cm,
      prefer_hrel = prefer_hrel,
      rsteiner = r_steiner,
      rin = r_in,
      rsteiner_se = r_steiner_se,
      rin_se = r_in_se,
      n = n,
      n_draws = n_draws,
      conf = conf,
      seed = seed,
      keep_draws = keep_draws,
      dsw_se = dsw_se,
      include_profile_uncertainty = include_profile_uncertainty,
      include_rmodel_residual = include_rmodel_residual
    ))
  }
  .SF_corr_value(profile = profile, geom = geom, models = models,
                 eval_mode = eval_mode, dsw = dsw,
                 active_depth = active_depth,
                 h_rel = h_rel, h_max = h_max, h_max_cm = h_max_cm,
                 prefer_hrel = prefer_hrel, n = n,
                 rsteiner = r_steiner, rin = r_in)
}

simulate_SF_ci <- function(
    fit_row,                 # one row: alpha, beta, (optional) vcov_log_11/12/22
    geom,
    models = NULL,           # needed only for SF_corr
    dsw_cm,
    dsw_se_cm = NULL,        # user-specified uncertainty in dsw (cm)
    method = c("circ","steiner","corr"),
    eval_mode = c("r","z"),
    n_draws = 2000,
    ci = c(0.025, 0.975),
    # optional: include uncertainty in convexity etc (later)
    seed = NULL
) {
  method <- match.arg(method)
  eval_mode <- match.arg(eval_mode)
  
  if (!is.null(seed)) set.seed(seed)
  
  # ---- draw alpha,beta ----
  alpha_hat <- fit_row$alpha
  beta_hat  <- fit_row$beta
  
  if (!is.finite(alpha_hat) || !is.finite(beta_hat) || alpha_hat <= 0 || beta_hat <= 0) {
    stop("fit_row must include positive alpha and beta.")
  }
  
  # default: if no vcov, do independent lognormal using se_alpha/se_beta if available,
  # else treat as fixed.
  have_vcov <- all(c("vcov_log_11","vcov_log_12","vcov_log_22") %in% names(fit_row)) &&
    all(is.finite(c(fit_row$vcov_log_11, fit_row$vcov_log_12, fit_row$vcov_log_22)))
  
  if (have_vcov) {
    V <- matrix(c(fit_row$vcov_log_11, fit_row$vcov_log_12,
                  fit_row$vcov_log_12, fit_row$vcov_log_22), 2, 2)
    # draw (log alpha, log beta) ~ MVN
    # (no external deps: use chol + rnorm)
    mu <- c(log(alpha_hat), log(beta_hat))
    L <- try(chol(V), silent = TRUE)
    if (inherits(L, "try-error")) {
      # fall back to diagonal
      V <- diag(pmax(diag(V), 0), 2)
      L <- chol(V + diag(1e-12, 2))
    }
    Z <- matrix(rnorm(n_draws * 2), ncol = 2)
    theta <- sweep(Z %*% t(L), 2, mu, "+")
    alpha_draw <- exp(theta[,1])
    beta_draw  <- exp(theta[,2])
    
  } else if (all(c("se_alpha","se_beta") %in% names(fit_row)) &&
             is.finite(fit_row$se_alpha) && is.finite(fit_row$se_beta) &&
             fit_row$se_alpha > 0 && fit_row$se_beta > 0) {
    
    # independent lognormal approximation
    sd_loga <- fit_row$se_alpha / alpha_hat
    sd_logb <- fit_row$se_beta  / beta_hat
    alpha_draw <- exp(rnorm(n_draws, mean = log(alpha_hat), sd = sd_loga))
    beta_draw  <- exp(rnorm(n_draws, mean = log(beta_hat),  sd = sd_logb))
    
  } else {
    alpha_draw <- rep(alpha_hat, n_draws)
    beta_draw  <- rep(beta_hat,  n_draws)
  }
  
  # ---- draw dsw ----
  # dsw must be >0; use truncated normal by rejection (simple + dependency-free)
  if (is.null(dsw_se_cm) || !is.finite(dsw_se_cm) || dsw_se_cm <= 0) {
    dsw_draw <- rep(dsw_cm, n_draws)
  } else {
    dsw_draw <- rnorm(n_draws, mean = dsw_cm, sd = dsw_se_cm)
    # simple truncation: resample negatives
    bad <- which(!is.finite(dsw_draw) | dsw_draw <= 0)
    it <- 0
    while (length(bad) > 0 && it < 50) {
      dsw_draw[bad] <- rnorm(length(bad), mean = dsw_cm, sd = dsw_se_cm)
      bad <- which(!is.finite(dsw_draw) | dsw_draw <= 0)
      it <- it + 1
    }
    # if still bad, clamp
    dsw_draw[dsw_draw <= 0] <- min(dsw_cm, 1e-6)
  }
  
  # ---- compute SF for each draw ----
  SF <- numeric(n_draws)
  
  for (i in seq_len(n_draws)) {
    prof <- make_custom_profile(
      mode = eval_mode,
      js_fun = function(x) dbeta(pmin(pmax(x, 0), 1), alpha_draw[i], beta_draw[i]),
      active_max = 1
    )
    
    if (method == "circ") {
      SF[i] <- SF_circ(profile = prof, geom = geom, eval_mode = eval_mode, dsw_cm = dsw_draw[i])
    } else if (method == "steiner") {
      SF[i] <- SF_steiner(profile = prof, geom = geom, eval_mode = eval_mode, dsw_cm = dsw_draw[i])
    } else {
      if (is.null(models)) stop("models is required for method='corr'")
      SF[i] <- SF_corr(profile = prof, geom = geom, models = models,
                       eval_mode = eval_mode, dsw_cm = dsw_draw[i])
    }
  }
  
  SF <- SF[is.finite(SF)]
  if (!length(SF)) stop("All SF draws were non-finite; check inputs.")
  
  list(
    SF_mean = mean(SF),
    SF_median = stats::median(SF),
    SF_ci = stats::quantile(SF, probs = ci, names = FALSE),
    draws = SF
  )
}


fit_beta_profiles <- function(
    data,
    basis = c("z","r"),
    level = c("tree","global"),          # pooled fits to compute (sensor always fit)
    plot = FALSE,
    
    # columns
    tree_col   = "tree",
    sensor_col = "sensor",
    depth_col  = "depth",
    js_col     = "sapflux",
    dsw_col    = "d_sw",                 # required for basis="z"
    Rcirc_col  = NULL,                   # required for basis="r" unless DBH_col given
    DBH_col    = NULL,                   # optional alternative to Rcirc_col for basis="r" (Rcirc=DBH/2)
    
    # fitting controls
    min_points = 2,                      # observed points per sensor (endpoint not counted)
    eps = 1e-4,
    n_grid = 200,
    endpoint_weight = 0.2,               # weight for the added (x≈1, Js=0) constraint
    amp_pool = c("mean","median"),       # how to pool sensor amplitude scalars into tree/global
    se = c("none","hessian"),
    verbose = FALSE
) {
  basis <- match.arg(basis)
  se <- match.arg(se)
  level <- unique(match.arg(level, several.ok = TRUE))
  amp_pool <- match.arg(amp_pool)
  
  .stopif <- function(cond, msg) if (!isTRUE(cond)) stop(msg, call. = FALSE)
  .warnif <- function(cond, msg) if (isTRUE(cond)) warning(msg, call. = FALSE)
  
  clamp01 <- function(x, eps = 1e-9) pmin(pmax(x, eps), 1 - eps)
  
  trapz <- function(x, y) {
    o <- order(x)
    x <- x[o]; y <- y[o]
    dx <- diff(x)
    sum((y[-1] + y[-length(y)]) * 0.5 * dx)
  }
  
  # ---- column checks ----
  need_cols <- c(tree_col, sensor_col, depth_col, js_col)
  if (basis == "z") need_cols <- c(need_cols, dsw_col)
  if (basis == "r") {
    if (is.null(Rcirc_col) && is.null(DBH_col)) {
      .stopif(FALSE, "basis='r' requires either Rcirc_col (depth/Rcirc) or DBH_col (Rcirc=DBH/2).")
    }
    if (!is.null(Rcirc_col)) need_cols <- c(need_cols, Rcirc_col)
    if (is.null(Rcirc_col) && !is.null(DBH_col)) need_cols <- c(need_cols, DBH_col)
  }
  missing_cols <- setdiff(need_cols, names(data))
  .stopif(length(missing_cols) == 0, paste("Missing required column(s):", paste(missing_cols, collapse = ", ")))
  
  df <- data
  df$.tree   <- df[[tree_col]]
  df$.sensor <- df[[sensor_col]]
  df$.depth  <- df[[depth_col]]
  df$.Js     <- df[[js_col]]
  
  .stopif(all(is.finite(df$.depth)), "All depth values must be finite.")
  .stopif(all(is.finite(df$.Js)), "All sapflux values must be finite (pre-filter NAs).")
  
  # ---- compute x in (0,1) ----
  if (basis == "z") {
    df$.dsw <- df[[dsw_col]]
    .stopif(all(is.finite(df$.dsw) & df$.dsw > 0), "All d_sw values must be finite and > 0.")
    
    # warn if dsw varies within (tree,sensor)
    dsw_by <- tapply(df$.dsw, interaction(df$.tree, df$.sensor, drop = TRUE),
                     function(v) length(unique(round(v, 6))))
    .warnif(any(dsw_by > 1),
            "Some (tree,sensor) groups have multiple d_sw values; using the first within each group.")
    
    df$.x <- df$.depth / df$.dsw
    bad <- df$.x > 1 + 1e-12
    if (any(bad)) {
      .warnif(TRUE, paste0("Dropping ", sum(bad), " rows with depth > d_sw (z>1)."))
      df <- df[!bad, , drop = FALSE]
    }
  } else {
    if (!is.null(Rcirc_col)) df$.Rcirc <- df[[Rcirc_col]] else df$.Rcirc <- df[[DBH_col]] / 2
    .stopif(all(is.finite(df$.Rcirc) & df$.Rcirc > 0), "All Rcirc values must be finite and > 0.")
    
    df$.x <- df$.depth / df$.Rcirc
    bad <- df$.x > 1 + 1e-12
    if (any(bad)) {
      .warnif(TRUE, paste0("Dropping ", sum(bad), " rows with depth > Rcirc (r>1)."))
      df <- df[!bad, , drop = FALSE]
    }
  }
  
  df$.x <- clamp01(df$.x, eps = eps)
  df$.key <- interaction(df$.tree, df$.sensor, drop = TRUE)
  
  # ---- fit beta PDF to observed intensity curve ----
  # Strategy:
  # 1) add endpoint x=1-eps, Js=0 (constraint)
  # 2) compute amplitude amp = ∫ Js dx (trapezoid)
  # 3) normalize: pdf_obs = Js / amp  (so ∫ pdf_obs dx = 1)
  # 4) fit dbeta(x; a,b) to pdf_obs by weighted SSE (weights ~ dx, plus endpoint_weight)
  fit_one_beta_pdf <- function(x, y, add_endpoint = TRUE) {
    ok <- is.finite(x) & is.finite(y)
    x <- x[ok]; y <- y[ok]
    y <- pmax(y, 0)
    
    n_obs <- length(x)
    if (n_obs < min_points) return(list(ok = FALSE, reason = "too_few_points"))
    
    w <- rep(1, length(x))
    
    if (add_endpoint) {
      x <- c(x, 1 - eps)
      y <- c(y, 0)
      w <- c(w, endpoint_weight)
    }
    
    # amplitude (area under raw Js vs x)
    amp <- trapz(x, y)
    if (!is.finite(amp) || amp <= 0) return(list(ok = FALSE, reason = "nonpositive_amp"))
    
    pdf_obs <- y / amp
    
    # weights: approximate equal influence per x-interval
    # compute dx-based weights for interior points
    o <- order(x)
    xs <- x[o]
    dx <- diff(xs)
    w_dx <- rep(0, length(xs))
    # assign each point half adjacent interval widths
    w_dx[1] <- dx[1] / 2
    w_dx[length(xs)] <- dx[length(dx)] / 2
    if (length(xs) > 2) {
      for (i in 2:(length(xs)-1)) w_dx[i] <- (dx[i-1] + dx[i]) / 2
    }
    # combine with provided endpoint weight (already in w)
    w_use <- w_dx
    # multiply by explicit weights vector (endpoint_weight already set)
    w_use <- w_use * w[o]
    # avoid zeros
    w_use[w_use == 0] <- min(w_use[w_use > 0], na.rm = TRUE)
    
    xs <- xs
    ys <- pdf_obs[o]
    
    # init guess using moments from pdf_obs weights
    mu0 <- sum(xs * ys * w_use) / sum(ys * w_use)
    mu0 <- min(max(mu0, 0.1), 0.9)
    K0 <- 8
    a0 <- mu0 * K0
    b0 <- (1 - mu0) * K0
    
    obj <- function(par) {
      loga <- par[1]; logb <- par[2]
      a <- exp(loga); b <- exp(logb)
      pred <- dbeta(xs, a, b)
      sum(w_use * (ys - pred)^2)
    }
    
    opt <- try(stats::optim(log(c(a0, b0)), obj, method = "BFGS",
                            hessian = (se == "hessian")), silent = TRUE)
    if (inherits(opt, "try-error") || opt$convergence != 0) {
      # fallback: smaller K init
      K1 <- 3
      a1 <- mu0 * K1
      b1 <- (1 - mu0) * K1
      opt <- try(stats::optim(log(c(a1, b1)), obj, method = "BFGS",
                              hessian = (se == "hessian")), silent = TRUE)
      if (inherits(opt, "try-error") || opt$convergence != 0) {
        return(list(ok = FALSE, reason = "optim_failed"))
      }
    }
    
    a <- exp(opt$par[1]); b <- exp(opt$par[2])
    mu <- a / (a + b)
    K  <- a + b
    
    out <- list(
      ok = TRUE,
      alpha = a,
      beta  = b,
      mu = mu,
      K  = K,
      amp = amp,          # amplitude scalar (area under raw Js curve)
      sse = opt$value
    )
    
    if (se == "hessian") {
      H <- opt$hessian
      if (!is.null(H) && all(is.finite(H))) {
        V <- try(solve(H), silent = TRUE)
        if (!inherits(V, "try-error") && all(diag(V) > 0)) {
          se_log <- sqrt(diag(V))
          out$se_alpha <- out$alpha * se_log[1]
          out$se_beta  <- out$beta  * se_log[2]
          # store covariance on log-scale
          out$vcov_log_11 <- V[1,1]
          out$vcov_log_12 <- V[1,2]
          out$vcov_log_22 <- V[2,2]
        }
      }
    }
    
    
    out
  }
  
  # ---- per-sensor fits ----
  keys <- levels(df$.key)
  sensor_fits <- vector("list", length(keys))
  
  for (i in seq_along(keys)) {
    k <- keys[i]
    sub <- df[df$.key == k, , drop = FALSE]
    f <- fit_one_beta_pdf(sub$.x, sub$.Js, add_endpoint = TRUE)
    
    sensor_fits[[i]] <- data.frame(
      tree   = sub$.tree[1],
      sensor = sub$.sensor[1],
      n_obs  = nrow(sub),
      basis  = basis,
      ok     = isTRUE(f$ok),
      reason = if (isTRUE(f$ok)) NA_character_ else f$reason,
      alpha  = if (isTRUE(f$ok)) f$alpha else NA_real_,
      beta   = if (isTRUE(f$ok)) f$beta  else NA_real_,
      mu     = if (isTRUE(f$ok)) f$mu else NA_real_,
      K      = if (isTRUE(f$ok)) f$K  else NA_real_,
      amp    = if (isTRUE(f$ok)) f$amp else NA_real_,   # amplitude scalar per sensor
      sse    = if (isTRUE(f$ok)) f$sse else NA_real_,
      se_alpha = if (!is.null(f$se_alpha)) f$se_alpha else NA_real_,
      se_beta  = if (!is.null(f$se_beta))  f$se_beta  else NA_real_,
      stringsAsFactors = FALSE
    )
  }
  sensor_fits <- do.call(rbind, sensor_fits)
  sensor_ok <- sensor_fits[sensor_fits$ok, , drop = FALSE]
  .stopif(nrow(sensor_ok) > 0, "No successful sensor fits; check data quality and min_points.")
  
  # grid for pooling
  x_grid <- seq(eps, 1 - eps, length.out = n_grid)
  pred_pdf <- function(a, b, x) dbeta(x, a, b)
  
  # helper to pool amplitude scalars
  pool_amp <- function(v) {
    v <- v[is.finite(v)]
    if (!length(v)) return(NA_real_)
    if (amp_pool == "median") stats::median(v) else mean(v)
  }
  
  # ---- per-tree pooled PDFs ----
  tree_fits <- NULL
  tree_curves <- NULL
  
  if ("tree" %in% level || "global" %in% level) {
    trees <- unique(sensor_ok$tree)
    tree_fits_list <- vector("list", length(trees))
    tree_curve_list <- vector("list", length(trees))
    
    for (i in seq_along(trees)) {
      tr <- trees[i]
      ss <- sensor_ok[sensor_ok$tree == tr, , drop = FALSE]
      if (nrow(ss) == 0) next
      
      # predicted PDFs on grid for each sensor
      mat <- vapply(seq_len(nrow(ss)), function(j) pred_pdf(ss$alpha[j], ss$beta[j], x_grid),
                    numeric(length(x_grid)))
      
      # average PDFs (sensor-equal). Average of PDFs integrates to 1 automatically
      pdf_mean <- if (is.matrix(mat)) rowMeans(mat) else mat
      
      # renormalize numerically to guard against tiny numeric drift
      pdf_mean <- pdf_mean / trapz(x_grid, pdf_mean)
      
      # refit a beta PDF to pooled pdf_mean
      ft <- {
        # treat pdf_mean as y and "amp=1" by construction:
        # We can reuse fit_one_beta_pdf by feeding y = pdf_mean and setting amp=1,
        # but fit_one_beta_pdf expects raw Js. Instead: fit alpha,beta directly on pdf_mean.
        xs <- x_grid
        ys <- pdf_mean
        # weights ~ dx
        dx <- diff(xs)
        w_use <- rep(0, length(xs))
        w_use[1] <- dx[1] / 2
        w_use[length(xs)] <- dx[length(dx)] / 2
        if (length(xs) > 2) for (k in 2:(length(xs)-1)) w_use[k] <- (dx[k-1] + dx[k]) / 2
        w_use[w_use == 0] <- min(w_use[w_use > 0])
        
        mu0 <- sum(xs * ys * w_use) / sum(ys * w_use)
        mu0 <- min(max(mu0, 0.1), 0.9)
        K0 <- 8
        a0 <- mu0 * K0
        b0 <- (1 - mu0) * K0
        
        obj <- function(par) {
          a <- exp(par[1]); b <- exp(par[2])
          pred <- dbeta(xs, a, b)
          sum(w_use * (ys - pred)^2)
        }
        
        opt <- try(stats::optim(log(c(a0, b0)), obj, method = "BFGS",
                                hessian = (se == "hessian")), silent = TRUE)
        if (inherits(opt, "try-error") || opt$convergence != 0) list(ok=FALSE, reason="optim_failed") else {
          a <- exp(opt$par[1]); b <- exp(opt$par[2])
          list(ok=TRUE, alpha=a, beta=b, mu=a/(a+b), K=a+b, sse=opt$value)
        }
      }
      
      tree_fits_list[[i]] <- data.frame(
        tree  = tr,
        n_sensors = nrow(ss),
        basis = basis,
        ok    = isTRUE(ft$ok),
        reason = if (isTRUE(ft$ok)) NA_character_ else ft$reason,
        alpha = if (isTRUE(ft$ok)) ft$alpha else NA_real_,
        beta  = if (isTRUE(ft$ok)) ft$beta  else NA_real_,
        mu    = if (isTRUE(ft$ok)) ft$mu else NA_real_,
        K     = if (isTRUE(ft$ok)) ft$K  else NA_real_,
        amp   = pool_amp(ss$amp),          # pooled amplitude scalar for that tree
        sse   = if (isTRUE(ft$ok)) ft$sse else NA_real_,
        stringsAsFactors = FALSE
      )
      
      tree_curve_list[[i]] <- data.frame(
        tree = tr,
        x = x_grid,
        pdf_mean = pdf_mean,
        pdf_fit  = if (isTRUE(ft$ok)) pred_pdf(ft$alpha, ft$beta, x_grid) else NA_real_
      )
    }
    
    tree_fits <- do.call(rbind, tree_fits_list)
    tree_curves <- do.call(rbind, tree_curve_list)
    tree_ok <- tree_fits[tree_fits$ok, , drop = FALSE]
    .stopif(nrow(tree_ok) > 0, "No successful tree fits; check sensor fits and pooling.")
  }
  
  # ---- global pooled PDF ----
  global_fit <- NULL
  global_curve <- NULL
  
  if ("global" %in% level) {
    tree_ok <- tree_fits[tree_fits$ok, , drop = FALSE]
    trees <- unique(tree_ok$tree)
    
    mat_tree <- vapply(seq_len(nrow(tree_ok)), function(j) pred_pdf(tree_ok$alpha[j], tree_ok$beta[j], x_grid),
                       numeric(length(x_grid)))
    pdf_global <- if (is.matrix(mat_tree)) rowMeans(mat_tree) else mat_tree
    pdf_global <- pdf_global / trapz(x_grid, pdf_global)
    
    # fit beta to global PDF
    xs <- x_grid; ys <- pdf_global
    dx <- diff(xs)
    w_use <- rep(0, length(xs))
    w_use[1] <- dx[1] / 2
    w_use[length(xs)] <- dx[length(dx)] / 2
    if (length(xs) > 2) for (k in 2:(length(xs)-1)) w_use[k] <- (dx[k-1] + dx[k]) / 2
    w_use[w_use == 0] <- min(w_use[w_use > 0])
    
    mu0 <- sum(xs * ys * w_use) / sum(ys * w_use)
    mu0 <- min(max(mu0, 0.1), 0.9)
    K0 <- 8
    a0 <- mu0 * K0
    b0 <- (1 - mu0) * K0
    
    obj <- function(par) {
      a <- exp(par[1]); b <- exp(par[2])
      pred <- dbeta(xs, a, b)
      sum(w_use * (ys - pred)^2)
    }
    opt <- try(stats::optim(log(c(a0, b0)), obj, method = "BFGS",
                            hessian = (se == "hessian")), silent = TRUE)
    
    if (!inherits(opt, "try-error") && opt$convergence == 0) {
      a <- exp(opt$par[1]); b <- exp(opt$par[2])
      global_fit <- data.frame(
        level = "global",
        n_trees = length(trees),
        basis = basis,
        ok = TRUE,
        alpha = a,
        beta  = b,
        mu = a/(a+b),
        K  = a+b,
        amp = pool_amp(tree_ok$amp),     # pooled amplitude scalar for global
        sse = opt$value,
        stringsAsFactors = FALSE
      )
      global_curve <- data.frame(
        x = x_grid,
        pdf_mean = pdf_global,
        pdf_fit  = pred_pdf(a, b, x_grid)
      )
    } else {
      global_fit <- data.frame(
        level="global", n_trees=length(trees), basis=basis, ok=FALSE,
        alpha=NA_real_, beta=NA_real_, mu=NA_real_, K=NA_real_, amp=NA_real_, sse=NA_real_,
        stringsAsFactors = FALSE
      )
      global_curve <- data.frame(x=x_grid, pdf_mean=pdf_global, pdf_fit=NA_real_)
    }
  }
  
  # ---- augmented output (same rows as input, plus pooled params) ----
  augmented <- data
  augmented$.tree <- data[[tree_col]]
  augmented$.sensor <- data[[sensor_col]]
  
  if (!is.null(tree_fits)) {
    tf <- tree_fits[, c("tree","alpha","beta","mu","K","amp","ok")]
    names(tf) <- c(".tree","alpha_tree","beta_tree","mu_tree","K_tree","amp_tree","tree_fit_ok")
    augmented <- merge(augmented, tf, by = ".tree", all.x = TRUE, sort = FALSE)
  }
  
  if (!is.null(global_fit) && isTRUE(global_fit$ok[1])) {
    augmented$alpha_global <- global_fit$alpha[1]
    augmented$beta_global  <- global_fit$beta[1]
    augmented$mu_global    <- global_fit$mu[1]
    augmented$K_global     <- global_fit$K[1]
    augmented$amp_global   <- global_fit$amp[1]
  }
  
  # ---- plotting ----
  if (isTRUE(plot)) {
    .stopif(requireNamespace("ggplot2", quietly = TRUE),
            "plot=TRUE requires ggplot2 installed.")
    
    # build per-row observed PDFs for plotting points:
    # pdf_obs = Js / amp_sensor
    amp_by_key <- setNames(sensor_ok$amp, interaction(sensor_ok$tree, sensor_ok$sensor, drop = TRUE))
    plot_df <- df
    plot_df$.amp_sensor <- amp_by_key[as.character(plot_df$.key)]
    plot_df$.pdf_obs <- plot_df$.Js / plot_df$.amp_sensor
    plot_df$tree <- plot_df$.tree
    
    xlab <- if (basis == "z") "z = depth / d_sw" else "r = depth / Rcirc"
    p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = .x, y = .pdf_obs, color = factor(tree))) +
      ggplot2::geom_point(alpha = 0.7) +
      ggplot2::labs(x = xlab, y = "Observed PDF (normalized within sensor)", color = "tree") +
      ggplot2::theme_bw()
    
    # per-tree fitted PDFs
    if (!is.null(tree_fits)) {
      tf_ok <- tree_fits[tree_fits$ok, , drop = FALSE]
      if (nrow(tf_ok) > 0) {
        tree_line_df <- do.call(rbind, lapply(seq_len(nrow(tf_ok)), function(i) {
          tr <- tf_ok$tree[i]
          data.frame(
            tree = tr,
            x = x_grid,
            y = pred_pdf(tf_ok$alpha[i], tf_ok$beta[i], x_grid)
          )
        }))
        p <- p + ggplot2::geom_line(data = tree_line_df,
                                    ggplot2::aes(x = x, y = y, color = factor(tree)),
                                    linewidth = 0.9, alpha = 0.9, inherit.aes = FALSE)
      }
    }
    
    # global curve
    if (!is.null(global_fit) && isTRUE(global_fit$ok[1])) {
      gdf <- data.frame(x = x_grid, y = global_curve$pdf_fit)
      p <- p + ggplot2::geom_line(data = gdf, ggplot2::aes(x = x, y = y),
                                  color = "black", linewidth = 1.6, inherit.aes = FALSE)
      
      lab <- sprintf("Global beta PDF: alpha=%.2f  beta=%.2f\nmu=%.2f  K=%.2f",
                     global_fit$alpha[1], global_fit$beta[1], global_fit$mu[1], global_fit$K[1])
      p <- p + ggplot2::annotate("label",
                                 x = 0.62,
                                 y = max(plot_df$.pdf_obs, na.rm = TRUE) * 0.95,
                                 label = lab, hjust = 0, size = 3.2)
    }
    
    print(p)
  }
  
  out <- list(
    basis = basis,
    sensor_fits = sensor_fits,     # includes amp per sensor
    tree_fits   = tree_fits,       # includes pooled amp per tree
    global_fit  = global_fit,      # includes pooled amp global
    augmented_data = augmented,
    curves = list(
      x_grid = x_grid,
      tree_curves = tree_curves,
      global_curve = global_curve
    )
  )
  
  if (verbose) {
    message("Sensor fits ok: ", sum(sensor_fits$ok), "/", nrow(sensor_fits))
    if (!is.null(tree_fits)) message("Tree fits ok: ", sum(tree_fits$ok), "/", nrow(tree_fits))
    if (!is.null(global_fit)) message("Global fit ok: ", isTRUE(global_fit$ok[1]))
  }
  
  out
}
