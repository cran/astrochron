\name{sedRamp}
\alias{sedRamp}
\title{Apply 'ramping' sedimentation rate model to convert time to stratigraphy}
\description{
Apply a linearly increasing (or decreasing) sedimentation rate model to convert time to stratigraphy.
}
\usage{
sedRamp(dat,srstart=0.01,srend=0.05,genplot=T,verbose=T)
}
\arguments{
  \item{dat}{Time series. First column should be time (in ka), second column should be data value.}
  \item{srstart}{Initial sedimentation rate (in m/ka).}
  \item{srend}{Final sedimentation rate (in m/ka).}
  \item{genplot}{Generate summary plots? (T or F)}
  \item{verbose}{Verbose output? (T or F)}
}

\value{
modeled stratigraphic series.
}

\examples{
\donttest{
# generate example series with 3 precession terms using function 'cycles'
# then convert from time to space using sedimentation rate that increases from 1 to 7 cm/ka
ex=sedRamp(cycles(),srstart=0.01,srend=0.07)
}
}