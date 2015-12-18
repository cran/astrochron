\name{timeOptSim}
\alias{timeOptSim}
\title{Monte Carlo Simulation for TimeOpt}
\description{
Perform Monte Carlo AR1 simulations to evaluate significance of TimeOpt results, as in Meyers (2015).
}
\usage{
timeOptSim(dat,sedrate=NULL,numsim=1000,rho=NULL,fit=1,flow=NULL,fhigh=NULL,
       roll=NULL,targetE=NULL,targetP=NULL,output=0,genplot=T,verbose=T)
}
\arguments{
  \item{dat}{Stratigraphic series for modulation assessment. First column should be depth or height (in meters), second column should be data value.}
  \item{sedrate}{Sedimentation rate for investigation (cm/ka).}
  \item{numsim}{Number of Monte Carlo AR1 simulations}
  \item{rho}{AR1 coefficient to use in simulations. By default this will be estimated from the stratigraphic series.}
  \item{fit}{Test for (1) precession amplitude modulations or (2) short eccentricity amplitude modulations? Option 2 is experimental!}
  \item{flow}{Low frequency cut-off for Taner bandpass (half power point; in cycles/ka)}
  \item{fhigh}{High frequency cut-off for Taner bandpass (half power point; in cycles/ka)}
  \item{roll}{Taner filter roll-off rate, in dB/octave.}
  \item{targetE}{A vector of eccentricity periods to evaluate (in ka). These must be in order of decreasing period, with a first value of 405 ka.}
  \item{targetP}{A vector of precession periods to evaluate (in ka). These must be in order of decreasing period.}
  \item{output}{Return: (0) nothing; (1) p-value; (2) output simulation r^2_opt results}
  \item{genplot}{Generate summary plots? (T or F)}
  \item{verbose}{Verbose output? (T or F)}
}

\seealso{ \code{\link{asm}}, \code{\link{eAsm}}, \code{\link{eAsmTrack}}, \code{\link{testPrecession}}, and \code{\link{timeOpt}} }

\references{
S.R. Meyers, 2015,
\emph{The evaluation of eccentricity-related amplitude modulation and bundling in paleoclimate data: An inverse approach for astrochronologic testing and time scale optimization}: Paleoceanography, 30, doi:10.1002/2015PA002850.
}

\examples{
\donttest{
# generate a test signal with precession and eccentricity
ex=etp(tmin=1,tmax=1000,dt=1,pWt=1,oWt=0,eWt=1,esinw=TRUE,genplot=FALSE,verbose=FALSE)
# convert to meters with sedimentation rate of 2 cm/kyr
ex[1]<-ex[1]*0.02
# evaluate precession modulations
timeOptSim(ex,sedrate=2,numsim=2000,fit=1)
}
}