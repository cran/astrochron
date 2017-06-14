\name{trackFreq}
\alias{trackFreq}
\title{Frequency-domain minimal tuning: Use interactive graphical interface and sorting to track frequency drift}
\description{
Frequency-domain minimal tuning: Use interactive graphical interface and sorting algorithm to track frequency drift.
}
\usage{
trackFreq(spec,threshold=NULL,pick=T,fmin=NULL,fmax=NULL,dmin=NULL,dmax=NULL,xmin=NULL,
          xmax=NULL,ymin=NULL,ymax=NULL,h=6,w=4,ydir=1,ncolors=100,genplot=T,verbose=T)
}
\arguments{
  \item{spec}{Time-frequency spectral results to evaluate. Must have the following format: column 
  1=frequency; remaining columns (2 to n)=power, amplitude or probability; titles for columns 2 to 
  n must be the location (depth or height). Note that this format is ouput by function eha.}
  \item{threshold}{Threshold level for filtering peaks. By default all peak maxima reported.}
  \item{pick}{Pick the peaks of interest using a graphical interface? (T or F). Only activated if genplot=T.}
  \item{fmin}{Minimum frequency for analysis.}
  \item{fmax}{Maximum frequency for analysis.}
  \item{dmin}{Minimum depth/height for analysis. NOT ACTIVATED YET!}
  \item{dmax}{Maximum depth/height for analysis. NOT ACTIVATED YET!}  
  \item{xmin}{Minimum frequency for PLOTTING.}
  \item{xmax}{Maximum frequency for PLOTTING.}
  \item{ymin}{Minimum depth/height for PLOTTING.}
  \item{ymax}{Maximum depth/height for PLOTTING.} 
  \item{h}{Height of plot in inches.}
  \item{w}{Width of plot in inches.}
  \item{ydir}{Direction for y-axis in plots (depth or height). -1 = values increase downwards (slower plotting!), 1 = values increase upwards.} 
  \item{ncolors}{Number of colors to use in plot.}
  \item{genplot}{Generate summary plots? (T or F)}
  \item{verbose}{Verbose output? (T or F)}
}

\seealso{ \code{\link{eha}} and \code{\link{traceFreq}} }

\examples{
\donttest{
# Check to see if this is an interactive R session, for compliance with CRAN standards.
# YOU CAN SKIP THE FOLLOWING LINE IF YOU ARE USING AN INTERACTIVE SESSION.
if(interactive()) {

# Generate example series with 3 terms using function 'cycles'.
# Then convert from time to space with sedimentation rate that increases from 1 to 5 cm/ka, using
# function 'sedramp'.
# Finally interpolate to median sampling interval using function 'linterp'.
dat=linterp(sedRamp(cycles(freqs=c(1/100,1/40,1/20),start=1,end=2500,dt=5)))

# EHA anlaysis, output probability results
out=eha(dat,output=4)

## Isolate peaks with probability >= 0.8
freq=trackFreq(out,0.8)

}
}
}