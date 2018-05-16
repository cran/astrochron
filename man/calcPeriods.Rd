\name{calcPeriods}
\alias{calcPeriods}
\title{Calculate eccentricity and precession periods in ka, given g and k in arcsec/yr}
\description{
Calculate eccentricity and precession periods in ka, given g and k in arcsec/yr.
}
\usage{
calcPeriods(g,k,output=1)
}
\arguments{
  \item{g}{Data frame or matrix with columns representing the fundamental frequencies: g1, g2, g3, g4, g5. Frequencies must be in arcsec/yr.}
  \item{k}{Data frame or vector with precession constant (frequency). Frequencies must be in arcsec/yr.}
  \item{output}{(1) return results as data frame, (2) return results as a numeric vector.}
}