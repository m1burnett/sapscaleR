test_that("make_trunk accepts diameter and radius", {
  by_diameter <- make_trunk(convexity = 0.9, diameter = 40)
  by_radius <- make_trunk(convexity = 0.9, radius = 20)

  expect_equal(by_diameter$radius, 20)
  expect_equal(by_radius$diameter, 40)
  expect_equal(by_diameter$convexity, by_radius$convexity)
})

test_that("beta curve can be applied on r, z, and d scales", {
  curve <- make_beta_curve(mu = 0.348, K = 3.66)
  trunk <- make_trunk(convexity = 0.85, diameter = 40)

  r_est <- SF_circ(curve, trunk, eval_mode = "r")
  z_est <- SF_circ(curve, trunk, eval_mode = "z", dsw = 10)
  d_est <- SF_circ(curve, trunk, eval_mode = "d", active_depth = 5)

  expect_true(is.finite(r_est))
  expect_true(is.finite(z_est))
  expect_true(is.finite(d_est))
  expect_gt(r_est, z_est)
  expect_gt(z_est, d_est)
})

test_that("SF_corr falls back without h_rel", {
  curve <- make_beta_curve(mu = 0.348, K = 3.66)
  trunk <- make_trunk(convexity = 0.8, diameter = 40)

  expect_warning(
    out <- SF_corr(curve, trunk, eval_mode = "z", dsw = 10),
    "h_rel not provided"
  )
  expect_true(is.numeric(out))
  expect_length(out, 1)
  expect_true(is.finite(out))
})

test_that("MC output has a consistent structure", {
  curve <- make_beta_curve(mu = 0.348, K = 3.66)
  trunk <- make_trunk(convexity = 0.85, diameter = 40)

  out <- SF_corr(
    curve,
    trunk,
    eval_mode = "z",
    dsw = 10,
    h_rel = 0.1,
    uncertainty = "mc",
    n_draws = 50,
    conf = 0.9,
    seed = 1
  )

  expect_named(
    out,
    c("estimate", "se", "ci", "conf", "uncertainty_sources",
      "n_draws", "n_valid", "n_invalid")
  )
  expect_equal(out$conf, 0.9)
  expect_true(all(c("r_steiner", "r_in") %in% out$uncertainty_sources))
  expect_true(out$ci[["lower"]] <= out$estimate)
  expect_true(out$estimate <= out$ci[["upper"]])
})

test_that("MC accepts uncertainty for user-supplied corrected radii", {
  curve <- make_beta_curve(mu = 0.348, K = 3.66)
  trunk <- make_trunk(convexity = 0.85, diameter = 40)

  out <- SF_corr(
    curve,
    trunk,
    eval_mode = "z",
    dsw = 10,
    r_steiner = 0.08,
    r_in = 0.4,
    r_steiner_se = 0.005,
    r_in_se = 0.01,
    uncertainty = "mc",
    n_draws = 50,
    seed = 2,
    keep_draws = TRUE
  )

  expect_true(is.finite(out$estimate))
  expect_equal(out$n_valid + out$n_invalid, out$n_draws)
  expect_true(all(c("r_steiner", "r_in") %in% out$uncertainty_sources))
  expect_gt(stats::sd(out$draws), 0)
})

test_that("custom profile uncertainty is included in MC draws and sources", {
  profile <- make_custom_profile(
    mode = "d",
    js_fun = function(d) pmax(0, 5 * exp(-d / 3)),
    js_se = 0.2,
    active_depth = 7
  )
  trunk <- make_trunk(convexity = 0.85, diameter = 40)

  out <- SF_circ(
    profile,
    trunk,
    eval_mode = "d",
    active_depth = 7,
    uncertainty = "mc",
    n_draws = 50,
    seed = 3,
    keep_draws = TRUE
  )

  expect_true("profile" %in% out$uncertainty_sources)
  expect_gt(stats::sd(out$draws), 0)
})

test_that("fit_sapflux_beta_profile carries profile uncertainty for MC", {
  path <- system.file("extdata", "avg.sfd.pigr.csv", package = "sapscaleR")
  skip_if(!nzchar(path), "example sap flux data not installed")
  avg_sfd <- read.csv(path)

  fit <- fit_sapflux_beta_profile(
    avg_sfd,
    scale = "z",
    tree_col = "Tree.id",
    sensor_col = "Sensor.id",
    depth_col = "T.depth.cm",
    flux_col = "Avg.SFD.cm.h",
    dbh_col = "DBH.cm",
    dsw_col = "swd"
  )

  expect_true(fit$method %in% c("nlme", "optim"))
  expect_true(is.matrix(fit$profile$vcov_log))
  expect_true(all(dim(fit$profile$vcov_log) == c(2, 2)))
  expect_true(is.finite(fit$profile$se_alpha))
  expect_true(is.finite(fit$profile$se_beta))
})
