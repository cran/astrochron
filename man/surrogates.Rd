\name{surrogates}
\alias{surrogates}
\title{Generate phase-randomized surrogate series as in Ebisuzaki (1997)}
\description{
Generate phase-randomized surrogate series as in Ebisuzaki (1997).
}
\usage{
surrogates(dat,nsim=1,preserveMean=T,std=T,genplot=T,verbose=T)
}
\arguments{
  \item{dat}{Data series with one or two columns. If two columns, first should be location (e.g., depth), second column should be data value.}
  \item{nsim}{Number of phase-randomized surrogate series to generate.}
  \item{preserveMean}{Should surrogate series have the same mean value as data series? (T or F)}
  \item{std}{Standardize results to guarantee equivalent variance as data series? (T or F)}
  \item{genplot}{Generate summary plots? Only applies if nsim=1. (T or F)}
  \item{verbose}{Verbose output? (T or F)}
}

\details{
This function will generate phase-randomized surrogate series as in Ebisuzaki (1997). It
is an R-translation of the Matlab code by V. Moron (see link below), with modifications and
additional features.
}

\references{
W. Ebisuzaki, 1997,
\emph{A Method to Estimate the Statistical Significance of a Correlation When the Data Are Serially Correlated}: Journal of Climate, v. 10, p. 2147-2153.

Matlab code by V. Moron: http://www.mathworks.com/matlabcentral/fileexchange/10881-weaclim/content/ebisuzaki.m

Original C-code by W. Ebisuzaki: http://www.ftp.cpc.ncep.noaa.gov/wd51we/random_phase/
}

\examples{
# generate example series with 3 precession terms and noise
ex <- cycles(start=0,end=500,noisevar=.0004,dt=5)

# generate phase-randomized surrogates 
ran_ex <- surrogates(ex,nsim=1)

# compare periodograms of data and surrogates
res1 <- periodogram(ex,padfac=0,output=1,genplot=FALSE)
res2 <- periodogram(ran_ex,padfac=0,output=1,genplot=FALSE)

pl(2)
plot(ex,type="l",main="black=original; red=surrogate")
lines(ran_ex,col="red",lty=4)
plot(res1[,1],res1[,2],type="l",lwd=2,main="black=original; red=surrogate",
     xlab="frequency",ylab="amplitude")
lines(res2[,1],res2[,2],col="red",lwd=2,lty=4)
}