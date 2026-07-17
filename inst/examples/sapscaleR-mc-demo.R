# Run from the repository root:
#   Rscript inst/examples/sapscaleR-mc-demo.R

if (requireNamespace("sapscaleR", quietly = TRUE)) {
  library(sapscaleR)
} else {
  source(file.path("R", "sapscaleR.R"))
}

trunk <- make_trunk(
  convexity = 0.9,
  diameter = 40
)

profile <- make_beta_curve(
  mu = 0.35,
  K = 8,
  se_alpha = 0.25,
  se_beta = 0.35
)

dsw <- 10

deterministic <- c(
  SF_circ = SF_circ(profile, trunk, eval_mode = "z", dsw = dsw),
  SF_steiner = SF_steiner(profile, trunk, eval_mode = "z", dsw = dsw),
  SF_corr_hrel = SF_corr(profile, trunk, eval_mode = "z", dsw = dsw, h_rel = 0.1),
  SF_corr_no_hrel = SF_corr(profile, trunk, eval_mode = "z", dsw = dsw),
  SF_corr_user_r = SF_corr(
    profile, trunk,
    eval_mode = "z",
    dsw = dsw,
    r_steiner = 0.25,
    r_in = 0.75
  )
)

print(round(deterministic, 2))

summarize_mc <- function(x) {
  c(
    estimate = x$estimate,
    se = x$se,
    lower = x$ci[["lower"]],
    upper = x$ci[["upper"]],
    n_valid = x$n_valid,
    n_invalid = x$n_invalid
  )
}

mc_corr_hrel <- SF_corr(
  profile, trunk,
  eval_mode = "z",
  dsw = dsw,
  h_rel = 0.1,
  uncertainty = "mc",
  n_draws = 500,
  seed = 1,
  dsw_se = 0.5
)

mc_corr_no_hrel <- SF_corr(
  profile, trunk,
  eval_mode = "z",
  dsw = dsw,
  uncertainty = "mc",
  n_draws = 500,
  seed = 1,
  dsw_se = 0.5
)

mc_circ_fixed_dsw <- SF_circ(
  profile, trunk,
  eval_mode = "z",
  dsw = dsw,
  uncertainty = "mc",
  n_draws = 500,
  seed = 2
)

mc_circ_uncertain_dsw <- SF_circ(
  profile, trunk,
  eval_mode = "z",
  dsw = dsw,
  uncertainty = "mc",
  n_draws = 500,
  seed = 2,
  dsw_se = 1
)

mc_summary <- rbind(
  corr_hrel = summarize_mc(mc_corr_hrel),
  corr_no_hrel = summarize_mc(mc_corr_no_hrel),
  circ_fixed_dsw = summarize_mc(mc_circ_fixed_dsw),
  circ_uncertain_dsw = summarize_mc(mc_circ_uncertain_dsw)
)

print(round(mc_summary, 2))
