\name{tune}
\alias{tune}
\title{Tune stratigraphic series}
\description{
Tune stratigraphic series from space to time, using specified control points
}
\usage{
tune(dat,controlPts,extrapolate=F,genplot=T,check=T,verbose=T)
}
\arguments{
  \item{dat}{Stratigraphic series for tuning. First column should be location (e.g., depth), second column should be data value.}
  \item{controlPts}{Tuning control points. A data frame or matrix containing two columns: depth, time}
  \item{extrapolate}{Extrapolate sedimentation rates above and below 'tuned' interval? (T or F)}
  \item{genplot}{Generate summary plots? (T or F)}
  \item{check}{Conduct compliance checks before processing? (T or F) In general this should be activated; the option is included for Monte Carlo simulation.}
  \item{verbose}{Verbose output? (T or F)}
}
\examples{
# generate example series with 3 precession terms using function 'cycles'
ex1=cycles()

# then convert from time to space using a sedimentation rate that increases from 1 to 7 cm/ka
ex2=sedRamp(ex1,srstart=0.01,srend=0.07)

# assemble tuning control points (this is the depth-time map)
controlPts=cbind(ex2[,1],ex1[,1])

# tune record
ex3=tune(ex2,controlPts=controlPts)
}