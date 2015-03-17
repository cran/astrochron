\name{traceFreq}
\alias{traceFreq}
\title{Frequency-domain minimal tuning: Use interactive graphical interface to trace frequency drift}
\description{
Frequency-domain minimal tuning: Use interactive graphical interface to trace frequency drift.
}
\usage{
traceFreq(spec,color=2,h=6,w=4,ydir=1,xmin=NULL,xmax=NULL,ymin=NULL,ymax=NULL,ncolors=100,
          pl=1,ln=F)
}
\arguments{
  \item{spec}{Time-frequency spectral results to evaluate. Must have the following format: column 
  1=frequency; remaining columns (2 to n)=power, amplitude or probability; titles for columns 2 to 
  n must be the location (depth or height). Note that this format is ouput by function eha.}
  \item{color}{Line color for tracing. 1 = transparent black; 2 = transparent white; 3 = transparent yellow}
  \item{h}{Height of plot in inches.}
  \item{w}{Width of plot in inches.}
  \item{ydir}{Direction for y-axis in plots (depth or height). -1 = values increase downwards (slower plotting!), 1 = values increase upwards.}
  \item{xmin}{Minimum spatial frequency to plot.}
  \item{xmax}{Maximum spatial frequency to plot.}
  \item{ymin}{Minimum depth/height to plot.}
  \item{ymax}{Maximum depth/height to plot.}  
  \item{ncolors}{Number of colors to use in plot.}
  \item{pl}{How do you want to represent the spatial frequency path?: 1=lines and points; 2=lines; 3=points}
  \item{ln}{Plot natural log of spectral results? (T or F)}
}

\seealso{ \code{\link{eha}} and \code{\link{trackFreq}} }

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

# EHA anlaysis, output amplitude results
out=eha(dat,output=3)

## Interactively track frequency drift
freq=traceFreq(out)

}
}
}