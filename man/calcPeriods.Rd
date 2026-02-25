\name{calcPeriods}
\alias{calcPeriods}
\title{Calculate eccentricity, obliquity and precession periods in ka, given g, s and k in arcsec/yr}
\description{
Calculate eccentricity, obliquity and precession periods in ka, given g, s and k in arcsec/yr.
}
\usage{
calcPeriods(g,s=NULL,k,opt=1,output=1)
}
\arguments{
  \item{g}{Required data frame or matrix with five columns representing the g-fundamental frequencies (g1, g2, g3, g4, g5). Frequencies must be in arcsec/yr.}
  \item{s}{Optional data frame or matrix with five or six columns representing the s-fundamental frequencies (five columns = s1, s2, s3, s4, s6; six columns = s1, s2, s3, s4, s5, s6). Frequencies must be in arcsec/yr. This is only required if you want to calculate obliquity periods.}
  \item{k}{Required data frame or vector with precession constant (frequency). Frequencies must be in arcsec/yr.}
  \item{opt}{(1) Calculate obliquity periods as in Laskar (1999, Table 6) (2) Calculate obliquity periods as in Malinverno and Meyers (2024, Table 1).}
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

\references{
Laskar, J., 1999,
\emph{The limits of Earth orbital calculations for geological time-scale use}: Philosophical Transactions, Mathematical, Physical and Engineering Sciences, 357, 1735-1759. 

Malinverno, A. and Meyers, S.R., 2024,
\emph{Bayesian estimation of past astronomical frequencies, lunar distance, and length of day from sediment cycles}: Geochemistry, Geophysics, Geosystems, 25, e2023GC011176.
}
