\name{calcPeriods}
\alias{calcPeriods}
\title{Calculate eccentricity, obliquity and precession periods in ka, given g, s and k in arcsec/yr}
\description{
Calculate eccentricity, obliquity and precession periods in ka, given g, s and k in arcsec/yr.
}
\usage{
calcPeriods(g,s=NULL,k,output=1)
}
\arguments{
  \item{g}{Required Data frame or matrix with columns representing the fundamental frequencies: g1, g2, g3, g4, g5. Frequencies must be in arcsec/yr.}
  \item{s}{Optional data frame or matrix with columns representing the fundamental frequencies: s1, s2, s3, s4, s5, s6. Frequencies must be in arcsec/yr. This is only required if you want to calculate obliquity periods.}
  \item{k}{Required data frame or vector with precession constant (frequency). Frequencies must be in arcsec/yr.}
  \item{output}{(1) return results as data frame, (2) return results as a numeric vector.}
}
\examples{
# calculate eccentricity and precession periods for one set of g's and k
gVal= c(5.579378,7.456665,17.366595,17.910194,4.257564)
kVal= 50.475838
calcPeriods(g=gVal, k=kVal)

# calculate eccentricity and precession periods for three sets of g's and k
gVal= matrix(c(5.579378,7.456665,17.366595,17.910194,4.257564,5.494302,7.452619,
      17.480760,18.348310,4.257451,5.531285,7.456848,17.320480,17.912240,4.257456),
      nrow=3,ncol=5,byrow=TRUE)

kVal= c(50.475838,51.280910,85.790450)

calcPeriods(g=gVal, k=kVal)
}