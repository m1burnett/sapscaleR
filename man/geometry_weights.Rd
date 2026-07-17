\name{g_circ}
\alias{g_circ}
\alias{g_steiner}
\alias{g_corr}
\title{Geometry weighting functions}
\usage{
g_circ(r)

g_steiner(r, convexity)

g_corr(r, convexity, rsteiner, rin, active_r_max = NULL)
}
\arguments{
\item{r}{Relative radial depth.}
\item{convexity}{Perimeter convexity.}
\item{rsteiner}{Relative depth where the Steiner correction stops being valid.}
\item{rin}{Relative inscribed circle radius.}
\item{active_r_max}{Legacy argument retained for compatibility.}
}
\description{
Low-level geometry weighting functions used by the sap flow estimators.
}
\value{Numeric vector of geometry weights.}
