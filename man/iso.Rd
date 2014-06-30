\name{iso}
\alias{iso}
\title{Isolate Data from a Specified Stratigraphic Interval}
\description{
Isolate a section of a uni- or multi-variate stratigraphic data set for further analysis
}
\usage{
iso(dat,xmin,xmax,col=2,genplot=T,verbose=T)
}
\arguments{
  \item{dat}{Data frame containing stratigraphic variable(s) of interest. First column must be location (e.g., depth).}
  \item{xmin}{Minimum depth/height/time for isolation. If xmin is not specified, it will be selected using a graphical interface.}
  \item{xmax}{Maximum depth/height/time for isolation. If xmax is not specified, it will be selected using a graphical interface.}
  \item{col}{If you are using the graphical interface to select xmin/xmax, which column would you like to plot? (default = 2).}
  \item{genplot}{Generate summary plots? (T or F)}
  \item{verbose}{Verbose output? (T or F)}
}

\seealso{ \code{\link{delPts}}, \code{\link{idPts}}, \code{\link{trim}} and \code{\link{trimAT}} }