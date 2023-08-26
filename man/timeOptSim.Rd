\name{timeOptSim}
\alias{timeOptSim}
\title{Monte Carlo simulation for TimeOpt}
\description{
Perform Monte Carlo AR1 simulations to evaluate significance of TimeOpt results, as in Meyers (2015).
}
\usage{
timeOptSim(dat,numsim=2000,rho=NULL,sedrate=NULL,sedmin=0.5,sedmax=5,numsed=100,
         linLog=1,limit=T,fit=1,fitModPwr=T,flow=NULL,fhigh=NULL,roll=NULL,
         targetE=NULL,targetP=NULL,detrend=T,ncores=2,output=0,genplot=T,
         check=T,verbose=T)
 }

\arguments{
  \item{dat}{Stratigraphic series for astrochronologic assessment. First column should be depth or height (in meters), second column should be data value.}
  \item{numsim}{Number of Monte Carlo AR1 simulations.}
  \item{rho}{AR1 coefficient to use in simulations. By default this will be estimated from the stratigraphic series.}
  \item{sedrate}{Sedimentation rate for investigation (cm/ka). This option is for compatibility with prior versions of timeOptSim. Please use sedmin, sedmax, numsed.}
  \item{sedmin}{Minimum sedimentation rate for investigation (cm/ka).}
  \item{sedmax}{Maximum sedimentation rate for investigation (cm/ka).}
  \item{numsed}{Number of sedimentation rates to investigate in optimization grid.}
  \item{linLog}{Use linear or logarithmic scaling for sedimentation rate grid spacing? (0=linear, 1=log)}
  \item{limit}{Limit evaluated sedimentation rates to region in which full target signal can be recovered? (T or F)}
  \item{fit}{Test for (1) precession amplitude modulation or (2) short eccentricity amplitude modulation?}
  \item{fitModPwr}{Include the modulation periods in the spectral fit? (T or F)}
  \item{flow}{Low frequency cut-off for Taner bandpass (half power point; in cycles/ka)}
  \item{fhigh}{High frequency cut-off for Taner bandpass (half power point; in cycles/ka)}
  \item{roll}{Taner filter roll-off rate, in dB/octave.}  
  \item{targetE}{A vector of eccentricity periods to evaluate (in ka). These must be in order of decreasing period, with a first value of 405 ka.}
  \item{targetP}{A vector of precession periods to evaluate (in ka). These must be in order of decreasing period.}
  \item{detrend}{Remove linear trend from data series? (T or F)}
  \item{ncores}{Number of cores to use for parallel processing. Must be >=2}
  \item{output}{Which results would you like to return to console? (0) no output; (1) p-value; (2) simulation r2 results}
  \item{genplot}{Generate summary plots? (T or F)}
  \item{check}{Conduct compliance checks before processing? (T or F). In general this should be activated.}
  \item{verbose}{Verbose output? (T or F)}
}

\details{
TimeOpt is an astronomical testing algorithm for untuned (spatial) stratigraphic data. 
The algorithm identifies the sedimentation rate(s) that simultaneously optimizes: 
(1) eccentricity amplitude modulations within the precession band, and (2) the 
concentration of spectral power at specified target astronomical periods.
The statistical significance of the r^2_opt is determined via Monte Carlo simulation using timeOptSim.

The present version of timeOptSim improves upon the original significance testing method of Meyers (2015),
by conducting simulations across the entire sedimentation grid. This approach more rigorously
protects against inflation of the p-value due to multiple testing. Parallel processing has been implemented
to address the greater computational demand that is required.

See timeOpt for more information on the basic methodology.
}

\seealso{ \code{\link{asm}}, \code{\link{eAsm}}, \code{\link{eAsmTrack}}, \code{\link{testPrecession}}, \code{\link{timeOpt}}, and \code{\link{timeOptPlot}} }

\references{
S.R. Meyers, 2015,
\emph{The evaluation of eccentricity-related amplitude modulations and bundling in paleoclimate data: An inverse approach for astrochronologic testing and time scale optimization}: Paleoceanography.
}

\examples{
 \dontrun{
# generate a test signal with precession and eccentricity
ex=etp(tmin=1,tmax=1000,dt=5,pWt=1,oWt=0,eWt=1,esinw=TRUE,genplot=FALSE,verbose=FALSE)
# convert to meters with sedimentation rate of 2 cm/kyr
ex[1]<-ex[1]*0.02
# evaluate with timeOptSim. be patient, this may take a while to run.
timeOptSim(ex,sedmin=0.5,sedmax=5,numsed=100)
 }
}