\name{predict_rsteiner}
\alias{predict_rsteiner}
\alias{predict_rin}
\alias{hrel_from_hmax}
\alias{sapscale_models}
\title{Predict corrected-model radii}
\usage{
predict_rsteiner(convexity, models = sapscale_models)

predict_rin(convexity, h_rel = NULL, models = sapscale_models,
  prefer_hrel = TRUE, warn_missing_hrel = TRUE)

hrel_from_hmax(h_max = NULL, diameter = NULL, h_max_cm = NULL,
  DBH_cm = NULL)

sapscale_models
}
\arguments{
\item{convexity}{Perimeter convexity.}
\item{models}{Model coefficient list.}
\item{h_rel}{Deepest concavity depth divided by diameter.}
\item{prefer_hrel}{Prefer the convexity-plus-\code{h_rel} model when available.}
\item{warn_missing_hrel}{Warn when falling back to the convexity-only model.}
\item{h_max, h_max_cm}{Absolute deepest concavity depth.}
\item{diameter, DBH_cm}{Trunk diameter.}
}
\description{
Predict \code{r_steiner} and \code{r_in} using the paper-derived beta-regression
models. If \code{h_rel} is missing, \code{predict_rin()} can fall back to the
convexity-only Eq. S3.1 model.
}
\value{Numeric predicted relative radius or model list.}
