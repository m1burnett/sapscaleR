\name{make_trunk}
\alias{make_trunk}
\alias{make_geom}
\title{Create a trunk geometry object}
\usage{
make_trunk(convexity = NULL, perimeter = NULL, diameter = NULL,
  radius = NULL, hull_perimeter = NULL, perimeter_cm = NULL,
  DBH_cm = NULL, Ph_cm = NULL, Rcirc_cm = NULL)
}
\arguments{
\item{convexity}{Perimeter convexity, usually hull perimeter divided by true perimeter.}
\item{perimeter}{True trunk perimeter.}
\item{diameter}{Circular-equivalent diameter.}
\item{radius}{Circular-equivalent radius.}
\item{hull_perimeter}{Convex hull perimeter.}
\item{perimeter_cm, DBH_cm, Ph_cm, Rcirc_cm}{Legacy unit-specific aliases.}
}
\description{
Creates the trunk object used by the sap flow estimators. Either \code{diameter}
or \code{radius} is required. \code{convexity} can be supplied directly or
computed from \code{perimeter} and \code{hull_perimeter}.
}
\value{A list containing convexity, diameter, radius, and legacy aliases.}
