\name{replEps}
\alias{replEps}
\title{Replace Values <= 0 with Smallest Positive Value}
\description{
Replace all variable values <= 0 with the smallest positive floating-point number (eps) that can be represented on machine. If first column is location ID (depth/height/time), it will not be processed.  Any number of variables (columns) permitted.
}
\usage{
replEps(dat,ID=T,genplot=T,verbose=T)
}
\arguments{
  \item{dat}{Data series to process. If location is included (e.g., depth), it should be in the first column.}
  \item{ID}{Is a location ID included in the first column? (T or F)}
  \item{genplot}{Generate summary plots? (T or F)}
  \item{verbose}{Verbose output? (T or F)}
}