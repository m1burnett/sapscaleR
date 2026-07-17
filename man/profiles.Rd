\name{make_beta_curve}
\alias{make_beta_curve}
\alias{make_beta_profile}
\alias{make_custom_profile}
\alias{beta_peak_muK}
\title{Create sap flux profile objects}
\usage{
make_beta_curve(mu, K, peak = 1, vcov_log = NULL,
  se_alpha = NULL, se_beta = NULL)

make_beta_profile(mode = NULL, mu, K, peak = 1, vcov_log = NULL,
  se_alpha = NULL, se_beta = NULL, active_max = NULL,
  active_depth = NULL)

make_custom_profile(mode = c("r", "z", "d"), js_fun,
  active_max = NULL, active_depth = NULL, js_se = NULL,
  js_se_fun = NULL, draw_min = 0)

beta_peak_muK(mu, K)
}
\arguments{
\item{mu}{Beta distribution mean on the 0-to-1 scale.}
\item{K}{Beta precision-like parameter where alpha = mu * K and beta = (1 - mu) * K.}
\item{peak}{Multiplier applied to the beta density.}
\item{vcov_log}{Optional 2 by 2 covariance matrix for log alpha and log beta.}
\item{se_alpha, se_beta}{Optional standard errors for alpha and beta.}
\item{mode}{Legacy or custom-profile scale: \code{"r"}, \code{"z"}, or \code{"d"}.}
\item{active_max, active_depth}{Finite depth reached by a \code{"d"} profile.}
\item{js_fun}{Function returning sap flux density at the supplied depth coordinate.}
\item{js_se}{Optional constant standard error for custom-profile MC draws.}
\item{js_se_fun}{Optional function returning standard errors at the supplied depth coordinate.}
\item{draw_min}{Optional lower bound applied to custom-profile MC draws. Defaults to \code{0}.}
}
\description{
\code{make_beta_curve()} creates a beta curve on the unit interval. The curve is
placed in the trunk later by \code{eval_mode} in \code{SF_circ()},
\code{SF_steiner()}, or \code{SF_corr()}. \code{make_beta_profile()} is a
legacy wrapper that stores the profile mode in the object. \code{make_custom_profile()}
wraps arbitrary user-provided profile functions. Custom profiles can optionally
carry simple normal uncertainty via \code{js_se} or \code{js_se_fun}; MC draws
sample sap flux density around \code{js_fun(x)} and apply \code{draw_min}.
}
\value{A profile object with a \code{js_fun} function.}
