\name{fit_sapflux_beta_profile}
\alias{fit_sapflux_beta_profile}
\title{Fit a beta sap flux curve}
\usage{
fit_sapflux_beta_profile(data, scale = c("z", "r", "d"),
  tree_col = "Tree.id", sensor_col = "Sensor.id",
  depth_col = "T.depth.cm", flux_col = "Avg.SFD.cm.h",
  dbh_col = "DBH.cm", dsw_col = NULL, r_col = NULL, z_col = NULL,
  active_depth = NULL, active_depth_cm = NULL,
  normalize = c("sensor_max", "none"), exclude_sensor_ids = NULL,
  add_endpoint = TRUE, endpoint_x = 0.999, endpoint_y = 0.001,
  clamp_y = c(0.001, 0.999), method = c("nlme", "optim"),
  eps = 1e-06)
}
\arguments{
\item{data}{Data frame containing sap flux observations.}
\item{scale}{Scale used to fit the curve: \code{"z"}, \code{"r"}, or \code{"d"}.}
\item{tree_col, sensor_col, depth_col, flux_col}{Column names.}
\item{dbh_col, dsw_col, r_col, z_col}{Optional columns used to compute scaled depth.}
\item{active_depth, active_depth_cm}{Finite active depth for \code{scale = "d"}.}
\item{normalize}{Whether to normalize each sensor by its maximum.}
\item{exclude_sensor_ids}{Optional sensor IDs to remove before fitting.}
\item{add_endpoint, endpoint_x, endpoint_y}{Endpoint constraint controls.}
\item{clamp_y}{Range used to clamp normalized sap flux density.}
\item{method}{Fit with \code{nlme} when available or fixed-effect \code{optim}.}
\item{eps}{Small value used to avoid exact beta-distribution endpoints.}
}
\description{
Fits a beta curve to radial sap flux density observations on \code{z}, \code{r},
or \code{d} scales. The returned \code{profile} is a beta curve on the unit
interval and is placed in trunks later by the SF estimator \code{eval_mode}.
The returned profile carries an approximate covariance matrix for
\code{log(alpha)} and \code{log(beta)} when this can be estimated from the
\code{nlme} fit or the fixed-effect optimizer. MC sap flow estimates use that
matrix to draw profile-parameter uncertainty.
}
\value{A list containing fitted parameters, optional uncertainty estimates, fit
data, and profile objects.}
