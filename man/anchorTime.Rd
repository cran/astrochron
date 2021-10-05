\name{anchorTime}
\alias{anchorTime}
\title{Anchor a floating astrochronology to a radioisotopic age}
\description{
Anchor a floating astrochronology to a radioisotopic age. The floating astrochronology is centered on a given ('floating') time datum and assigned the 'anchored' age.
}
\usage{
anchorTime(dat,time,age,timeDir=1,flipOut=F,genplot=T,verbose=T)
}
\arguments{
  \item{dat}{Stratigraphic series. First column should be floating time scale, second column should be data value.}
  \item{time}{'Floating' time datum to center record on. Units should be ka.}
  \item{age}{Radioisotopic age (or othwerwise) for anchoring at floating 'time' datum. Units should be ka.}
  \item{timeDir}{Direction of 'floating' time in input record; 1 = elapsed time towards present; 2 = elapsed time away from present}
  \item{flipOut}{Flip the output (sort so the ages are presented in decreasing order)? (T or F)}
  \item{genplot}{Generate summary plots? (T or F)}
  \item{verbose}{Verbose output? (T or F)}
}
