\name{repl0}
\alias{repl0}
\title{Replace values < 0 with 0}
\description{
Replace all variable values < 0 with 0. If first column is location ID (depth/height/time), it will not be processed.  Any number of variables (columns) permitted.
 
}
\usage{
repl0(dat,ID=T,genplot=T,verbose=T)
}
\arguments{
  \item{dat}{Data series to process. If location is included (e.g., depth), it should be in the first column.}
  \item{ID}{Is a location ID included in the first column? (T or F)}
  \item{genplot}{Generate summary plots? (T or F)}
  \item{verbose}{Verbose output? (T or F)}
}