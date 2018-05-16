\name{timeOptPlot}
\alias{timeOptPlot}
\title{TimeOptPlot: Generate summary figure for TimeOpt analyses}
\description{
TimeOptPlot: Generate summary figure for TimeOpt analyses.
}
\usage{
timeOptPlot(dat=NULL,res1=NULL,res2=NULL,simres=NULL,fit=1,fitModPwr,flow=NULL,
            fhigh=NULL,roll=NULL,targetE=NULL,targetP=NULL,xlab="Depth (m)",
            ylab="Proxy Value",fitR=NULL,verbose=T)
}
\arguments{
  \item{dat}{Stratigraphic series used for astrochronologic assessment. First column should be depth or height (in meters), second column should be data value.}
  \item{res1}{Data frame containing TimeOpt results: sedimentation rate grid, r^2_envelope, r^2_power, r^2_opt.}
  \item{res2}{Data frame containing the optimal-fitted time series, bandpassed series, envelope, and reconstructed eccentricity model.}
  \item{simres}{Data frame containing the r^2_opt value for each Monte Carlo simulation.}
  \item{fit}{Test for (1) precession amplitude modulation or (2) short eccentricity amplitude modulation?}
  \item{fitModPwr}{Include the modulation periods in the spectral fit? (T or F)}
  \item{flow}{Low frequency cut-off for Taner bandpass (half power point; in cycles/ka).}
  \item{fhigh}{High frequency cut-off for Taner bandpass (half power point; in cycles/ka).}
  \item{roll}{Taner filter roll-off rate, in dB/octave.}
  \item{targetE}{A vector of eccentricity periods to evaluate (in ka). These must be in order of decreasing period, with a first value of 405 ka.}
  \item{targetP}{A vector of precession periods to evaluate (in ka). These must be in order of decreasing period.}
  \item{xlab}{Label for the depth/height axis.}
  \item{ylab}{Label for proxy variable evaluated.}
  \item{fitR}{The r2_opt value at the optimal sedimentation rate.}
  \item{verbose}{Verbose output? (T or F)}
}

\seealso{ \code{\link{asm}}, \code{\link{eAsmTrack}}, \code{\link{testPrecession}}, \code{\link{timeOpt}}, and \code{\link{timeOptSim}} }

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
res1=timeOpt(ex,sedmin=0.5,sedmax=5,numsed=100,fit=1,output=1)
res2=timeOpt(ex,sedmin=0.5,sedmax=5,numsed=100,fit=1,output=2)
simres=timeOptSim(ex,sedrate=2,numsim=2000,fit=1,output=2)
timeOptPlot(ex,res1,res2,simres,flow=0.035,fhigh=0.065,roll=10^3,
 targetE=c(405.6795,130.719,123.839,98.86307,94.87666),
 targetP=c(23.62069,22.31868,19.06768,18.91979),xlab="Depth (m)",
 ylab="Value",fitR=0.832,verbose=T)
}
}