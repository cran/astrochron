\name{makeNoise}
\alias{makeNoise}
\title{Generate noise surrogates from a theoretical power spectrum}
\description{
Generate noise surrogates from a theoretical power spectrum.
}
\usage{
makeNoise(S,dt=1,mean=0,sdev=1,addPt=F,nsim=1,genplot=T,verbose=T)

}
\arguments{
  \item{S}{Vector or 1-D data frame containing the theoretical power, from f(0) to the Nyquist frequency.}
  \item{dt}{Sampling interval for surrogate series}
  \item{mean}{Mean value for surrogate series}
  \item{sdev}{Standard deviation for surrogate series}
  \item{addPt}{Did you add a Nyquist frequency? (T or F)}  
  \item{nsim}{Number of surrogate series to generate}
  \item{genplot}{generate summary plots (T or F)}
  \item{verbose}{verbose output (T or F)}
}
\details{
These simulations use the random number generator of Matsumoto and Nishimura [1998]. 
The algorithm of Timmer and Konig (1995) is employed to generate surrogates from
any arbitrary theoretical power spectrum. See examples section below.
}
\references{
M. Matsumoto, and T. Nishimura, (1998), 
\emph{Mersenne Twister: A 623-dimensionally equidistributed uniform pseudo-random number generator}, ACM Transactions on Modeling and Computer Simulation, 8, 3-30.

J. Timmer and K. Konig (1995), 
\emph{On Generating Power Law Noise}, Astronomy and Astrophysics: v. 300, p. 707-710.
}

\examples{
# create theoretical AR1 spectrum, using rho of 0.8
rho=0.8
freq=seq(0,.5,by=0.005)
Nyq=max(freq)
AR1 = (1-(rho^2)) / (  1 - (2*rho*cos(pi*freq/Nyq)) + (rho^2) )
plot(freq,AR1,type="l")

# make noise surrogates from the theoretical AR1 spectrum
makeNoise(AR1)
}