\name{pwrLaw}
\alias{pwrLaw}
\title{Generate power law (1/f) noise surrogates}
\description{
Generate power law (1/f) noise surrogates, following the algorithm of Timmer and Konig (1995).
}
\usage{
pwrLaw(npts=1024,dt=1,mean=0,sdev=1,beta=2,fcut=0,nsim=1,genplot=T,verbose=T)

}
\arguments{
  \item{npts}{number of data points for 1/f surrogate time series}
  \item{dt}{sampling interval}
  \item{mean}{mean value for 1/f surrogate series}
  \item{sdev}{standard deviation for 1/f surrogate series}
  \item{beta}{power law coefficient. Positive number will yield a negative slope.}
  \item{fcut}{frequency cutoff: below this frequency a plateau will be modeled. Set to zero (default) for no plateau.}
  \item{nsim}{Number of surrogate series to generate}
  \item{genplot}{generate summary plots (T or F)}
  \item{verbose}{verbose output (T or F)}
}
\details{
These simulations use the random number generator of Matsumoto and Nishimura (1998). 
Power law noise series are generated following the algorithm of Timmer and Konig (1995).

}
\references{
M. Matsumoto, and T. Nishimura, (1998), 
\emph{Mersenne Twister: A 623-dimensionally equidistributed uniform pseudo-random number generator}, ACM Transactions on Modeling and Computer Simulation, 8, 3-30.

J. Timmer and K. Konig (1995), 
\emph{On Generating Power Law Noise}, Astronomy and Astrophysics: v. 300, p. 707-710.
}