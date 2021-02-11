\name{eTimeOpt}
\alias{eTimeOpt}
\title{eTimeOpt: Evolutive implementation of TimeOpt (Meyers, 2015; Meyers, 2019)}
\description{
eTimeOpt: Evolutive implementation of TimeOpt (Meyers, 2015; Meyers, 2019).
}
\usage{
eTimeOpt(dat,win=dt*100,step=dt*10,sedmin=0.5,sedmax=5,numsed=100,linLog=1,
    limit=T,fit=1,fitModPwr=T,flow=NULL,fhigh=NULL,roll=NULL,targetE=NULL,targetP=NULL,
    detrend=T,ydir=1,output=1,genplot=T,check=T,verbose=1)
}
\arguments{
  \item{dat}{Stratigraphic series for astrochronologic assessment. First column should be depth or height (in meters), second column should be data value.}
  \item{win}{Window size, in meters.}
  \item{step}{Step size for moving window, in meters.}
  \item{sedmin}{Minimum sedimentation rate for investigation (cm/ka).}
  \item{sedmax}{Maximum sedimentation rate for investigation (cm/ka).}
  \item{numsed}{Number of sedimentation rates to investigate in optimization grid.}
  \item{linLog}{Use linear or logarithmic scaling for sedimentation rate grid spacing? (0=linear, 1=log; default value  is 1)}
  \item{limit}{Limit evaluated sedimentation rates to region in which full target signal can be recovered? (T or F).}
  \item{fit}{Test for (1) precession amplitude modulation or (2) short eccentricity amplitude modulation?}
  \item{fitModPwr}{Include the modulation periods in the spectral fit? (T or F)}
  \item{flow}{Low frequency cut-off for Taner bandpass (half power point; in cycles/ka)}
  \item{fhigh}{High frequency cut-off for Taner bandpass (half power point; in cycles/ka)}
  \item{roll}{Taner filter roll-off rate, in dB/octave.}  
  \item{targetE}{A vector of eccentricity periods to evaluate (in ka). These must be in order of decreasing period, with a first value of 405 ka.}
  \item{targetP}{A vector of precession periods to evaluate (in ka). These must be in order of decreasing period.}
  \item{detrend}{Remove linear trend from data series? (T or F)}
  \item{ydir}{Direction for y-axis in plots (depth,height,time). -1 = values increase downwards (slower plotting), 1 = values increase upwards}
  \item{output}{Which results would you like to return to the console? (0) no output; (1) everything, (2) r^2_envelope, (3) r^2_power, (4) r^2_opt}
  \item{genplot}{Generate summary plots? (T or F)}
  \item{check}{Conduct compliance checks before processing? (T or F) In general this should be activated; the option is included for Monte Carlo simulation.}
  \item{verbose}{Verbose output? (0=nothing, 1=minimal, 2=a little more, 3=everything!)}
}

\references{
S.R. Meyers, 2015,
\emph{The evaluation of eccentricity-related amplitude modulations and bundling in 
paleoclimate data: An inverse approach for astrochronologic testing and time scale 
optimization}: Paleoceanography, v.30, 1625-1640.

S.R. Meyers, 2019, 
\emph{Cyclostratigraphy and the problem of astrochronologic testing}: Earth-Science Reviews
v.190, 190-223.
}

\seealso{ \code{\link{tracePeak}},\code{\link{trackPeak}},\code{\link{timeOpt}},\code{\link{timeOptSim}}, and \code{\link{eTimeOptTrack}}}

\examples{
\dontrun{
# generate a test signal with precession and eccentricity
ex=cycles(freqs=c(1/405.6795,1/130.719,1/123.839,1/98.86307,1/94.87666,1/23.62069, 
1/22.31868,1/19.06768,1/18.91979),end=4000,dt=5)

# convert to meters with a linearly increasing sedimentation rate from 0.01 m/kyr to 0.03 m/kyr
ex=sedRamp(ex,srstart=0.01,srend=0.03)

# interpolate to median sampling interval
ex=linterp(ex)

# evaluate precession & eccentricity power, and precession modulations
res=eTimeOpt(ex,win=20,step=1,fit=1,output=1)

# extract the optimal fits for the power optimization
sedrates=eTimeOptTrack(res[2])

# extract the optimal fits for the envelope*power optimization
sedrates=eTimeOptTrack(res[3])

# you can also interactively track the results using functions 'trackPeak' and 'tracePeak'
#  evaluate the results from the power optimization
sedrates=tracePeak(res[2])
sedrates=trackPeak(res[2])

#  evaluate the results from the envelope*power optimization
sedrates=tracePeak(res[3])
sedrates=trackPeak(res[3])

# evaluate precession & eccentricity power, and short-eccentricity modulations
eTimeOpt(ex,win=20,step=1,fit=2,output=0)
}
}