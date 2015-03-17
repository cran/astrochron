\name{linage}
\alias{linage}
\title{Tune stratigraphic series to an astronomical target using graphical interface}
\description{
Tune stratigraphic series to an astronomical target using graphical interface similar to Analyseries 'Linage' routine (Paillard et al, 1996).
}
\usage{
linage(dat,target,extrapolate=F,xmin=NULL,xmax=NULL,tmin=NULL,tmax=NULL,plotype=1,
       output=1,genplot=T)
}
\arguments{
  \item{dat}{Stratigraphic series for tuning, with two columns. First column is depth/height.}
  \item{target}{Astronomical tuning target series. First column is time.}
  \item{extrapolate}{Extrapolate sedimentation rates above and below 'tuned' interval? (T or F)}
  \item{xmin}{Minimum height/depth to plot.}
  \item{xmax}{Maximum height/depth to plot.}
  \item{tmin}{Minimum time value to plot.}
  \item{tmax}{Maximum time value to plot.}  
  \item{plotype}{Type of plot to generate: 1= points and lines, 2 = points, 3 = lines}
  \item{output}{Return which of the following? 1 = tuned stratigraphic series; 2 = age control points; 3 = tuned stratigraphic series and age control points}
  \item{genplot}{Generate additional summary plots (tuned record, time-space map, sedimentation rates)? (T or F)}

}

\references{
Paillard, D., L. Labeyrie and P. Yiou, 1996),
\emph{Macintosh program performs time-series analysis}: Eos Trans. AGU, v. 77, p. 379.
}

\examples{
\donttest{
# Check to see if this is an interactive R session, for compliance with CRAN standards.
# YOU CAN SKIP THE FOLLOWING LINE IF YOU ARE USING AN INTERACTIVE SESSION.
if(interactive()) {

# generate example series with 3 precession terms and noise using function 'cycles'
# then convert from time to space using sedimentation rate that increases from 1 to 7 cm/ka
ex=sedRamp(cycles(start=1,end=400, dt=2,noisevar=.00005),srstart=0.01,srend=0.07)

# create astronomical target series
targ=cycles(start=1,end=400,dt=2)

## manually tune
tuned=linage(ex,targ)

## should you need to flip the direction of the astronomical target series, use function 'cb':
tuned=linage(ex,cb(targ[1]*-1,targ[2]))

}

}
}