\name{SF_circ}
\alias{SF_circ}
\alias{SF_steiner}
\alias{SF_corr}
\title{Whole-trunk sap flow estimators}
\usage{
SF_circ(profile, geom, eval_mode = c("r", "z", "d"), dsw = NULL,
  dsw_cm = NULL, active_depth = NULL, n = 2001,
  uncertainty = c("none", "mc"), n_draws = 2000, conf = 0.95,
  seed = NULL, keep_draws = FALSE, dsw_se = NULL, dsw_se_cm = NULL,
  include_profile_uncertainty = TRUE)

SF_steiner(profile, geom, eval_mode = c("r", "z", "d"), dsw = NULL,
  dsw_cm = NULL, active_depth = NULL, n = 2001,
  uncertainty = c("none", "mc"), n_draws = 2000, conf = 0.95,
  seed = NULL, keep_draws = FALSE, dsw_se = NULL, dsw_se_cm = NULL,
  include_profile_uncertainty = TRUE)

SF_corr(profile, geom, models = sapscale_models, eval_mode = c("r", "z", "d"),
  dsw = NULL, dsw_cm = NULL, active_depth = NULL, h_rel = NULL,
  h_max = NULL, h_max_cm = NULL, r_steiner = NULL, r_in = NULL,
  r_steiner_se = NULL, r_in_se = NULL, prefer_hrel = TRUE, n = 2001,
  uncertainty = c("none", "mc"), n_draws = 2000, conf = 0.95,
  seed = NULL, keep_draws = FALSE, dsw_se = NULL, dsw_se_cm = NULL,
  include_profile_uncertainty = TRUE, include_rmodel_residual = TRUE)
}
\arguments{
\item{profile}{Profile object from \code{make_beta_curve()} or \code{make_custom_profile()}.}
\item{geom}{Trunk object from \code{make_trunk()}.}
\item{eval_mode}{How to place a beta curve in the trunk: relative radius \code{"r"}, relative sapwood depth \code{"z"}, or absolute depth \code{"d"}.}
\item{dsw, dsw_cm}{Sapwood depth for \code{eval_mode = "z"}.}
\item{active_depth}{Depth reached by a beta curve for \code{eval_mode = "d"}.}
\item{n}{Number of radial integration grid points.}
\item{uncertainty}{Return a numeric estimate (\code{"none"}) or Monte Carlo summary (\code{"mc"}).}
\item{n_draws}{Number of Monte Carlo draws when \code{uncertainty = "mc"}.}
\item{conf}{Confidence level for the Monte Carlo quantile interval. Defaults to \code{0.95}.}
\item{seed}{Optional random seed for Monte Carlo draws.}
\item{keep_draws}{Whether to include the raw valid draws in the returned Monte Carlo list.}
\item{dsw_se, dsw_se_cm}{Standard error of sapwood depth for Monte Carlo draws.}
\item{include_profile_uncertainty}{Whether to draw beta profile parameters when uncertainty information is present.}
\item{models}{Model coefficient list, usually \code{sapscale_models}.}
\item{h_rel, h_max, h_max_cm}{Relative or absolute deepest concavity depth for the \code{r_in} model.}
\item{r_steiner, r_in}{Optional user-supplied corrected-model radii. These
can be measured values or values computed from the user's own functions before
calling \code{SF_corr()}. Supplied values are treated as fixed inputs in Monte
Carlo mode unless \code{r_steiner_se} or \code{r_in_se} is supplied.}
\item{r_steiner_se, r_in_se}{Optional standard errors for user-supplied
\code{r_steiner} and \code{r_in}. When \code{uncertainty = "mc"}, supplied
radii with positive SEs are drawn from normal distributions and invalid
geometries are resampled.}
\item{prefer_hrel}{Prefer the convexity-plus-\code{h_rel} \code{r_in} model when possible.}
\item{include_rmodel_residual}{Whether MC draws include beta-regression residual scatter for \code{r_steiner} and \code{r_in}.}
}
\description{
Estimate whole-trunk sap flow under the circular model, the shallow-depth
Steiner correction, or the corrected model for irregular trunks.
}
\value{
If \code{uncertainty = "none"}, a numeric sap flow estimate. If
\code{uncertainty = "mc"}, a list with estimate, standard error, confidence
interval, \code{uncertainty_sources}, draw counts, and optionally the raw
draws. \code{uncertainty_sources} reports which uncertainty components were
active, such as \code{"profile"}, \code{"dsw"}, \code{"r_steiner"}, and
\code{"r_in"}.
}
