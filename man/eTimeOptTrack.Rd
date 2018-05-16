\name{eTimeOptTrack}
\alias{eTimeOptTrack}
\title{Track eTimeOpt r2 maxima}
\description{
Track eTimeOpt r2 maxima.
}
\usage{
eTimeOptTrack(res,threshold=0,ydir=-1,genplot=T,verbose=T)
}
\arguments{
  \item{res}{eTimeOpt r2 results. Must have the following format: column 
  1=sedimentation rate; remaining columns (2 to n)=r2; titles for columns 2 to 
  n must be the location (depth or height). Note that this format is ouput by function eTimeOpt.}
  \item{threshold}{Threshold r2-value for analysis and plotting.}
  \item{ydir}{Direction for y-axis in plots (depth or height). -1 = values increase downwards (slower plotting!), 1 = values increase upwards.}
  \item{genplot}{Generate summary plots? (T or F)}
  \item{verbose}{Verbose output? (T or F)}
}

\details{
Please see function eTimeOpt for details.
}

\seealso{ \code{\link{timeOpt}}, and \code{\link{eTimeOpt}} }

\examples{
\donttest{
# generate a test signal with precession and eccentricity
ex=cycles(freqs=c(1/405.6795,1/130.719,1/123.839,1/98.86307,1/94.87666,1/23.62069, 
1/22.31868,1/19.06768,1/18.91979),end=4000,dt=5)

# convert to meters with a linearly increasing sedimentation rate from 0.01 m/kyr to 0.03 m/kyr
ex=sedRamp(ex,srstart=0.01,srend=0.03)

# interpolate to median sampling interval
ex=linterp(ex)

# evaluate precession & eccentricity power, and precession modulations
res=eTimeOpt(ex,win=20,step=.5,fit=1,output=1)

# extract the optimal fits for the power optimization
sedrates=eTimeOptTrack(res[2])

# extract the optimal fits for the envelope*power optimization
sedrates=eTimeOptTrack(res[3])

# you can also interactively track the results using functions 'trackPeak' and 'tracePeak'
#  evaluate the results from the power optimization
sedrates=trackPeak(res[2])
sedrates=tracePeak(res[2])

#  evaluate the results from the envelope*power optimization optimization
sedrates=trackPeak(res[3])
sedrates=tracePeak(res[3])
}
}