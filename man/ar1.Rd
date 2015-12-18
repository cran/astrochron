\name{ar1}
\alias{ar1}
\title{Generate AR(1) surrogates}
\description{
Generate AR(1) surrogates. Implement shuffling algorithm of Meyers (2012) if desired.
}
\usage{
ar1(npts=1024,dt=1,mean=0,sdev=1,rho=0.9,shuffle=F,nsim=1,genplot=T,verbose=T)

}
\arguments{
  \item{npts}{number of time series data points}
  \item{dt}{sampling interval}
  \item{mean}{mean value for AR1 surrogate series}
  \item{sdev}{standard deviation for AR1 surrogate series}
  \item{rho}{AR(1) coefficient}
  \item{shuffle}{Apply secondary shuffle of Gaussian deviates before AR modeling}
  \item{nsim}{Number of AR1 surrogate series to generate}
  \item{genplot}{generate summary plots (T or F)}
  \item{verbose}{verbose output (T or F)}
}
\details{
These simulations use the random number generator of Matsumoto and Nishimura [1998]. 
If shuffle = T, the algorithm from Meyers (2012, pg. 11) is applied: (1) two
sets of random sequences of same the length are generated, (2) the first random 
sequence is then sorted, and finally (3) the permutation vector of the sorted sequence 
is used to reorder the second random number sequence. This is done to guard against
potential shortcomings in random number generation that are specific to spectral estimation.
}
\references{
S.R. Meyers, 2012,
\emph{Seeing red in cyclic stratigraphy: Spectral noise estimation for astrochronology}: Paleoceanography, v. 27, PA3328.
}