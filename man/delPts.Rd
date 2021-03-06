\name{delPts}
\alias{delPts}
\title{Interactively delete points in plot}
\description{
Interactively delete points in x,y plot.
}
\usage{
delPts(dat,del=NULL,cols=c(1,2),ptsize=1,xmin=NULL,xmax=NULL,ymin=NULL,ymax=NULL,
       plotype=1,genplot=T,verbose=T)
}
\arguments{
  \item{dat}{Data frame containing stratigraphic variable(s) of interest. Any number of columns permitted.}
  \item{del}{A vector of indices indicating points to delete. If specified, the interactive plot is disabled.}
  \item{cols}{If you are using the graphical interface, which columns would you like to plot? (default = 1 & 2).}
  \item{ptsize}{Size of plotted points.}
  \item{xmin}{Minimum x-value (column 1) to plot}
  \item{xmax}{Maximum x-value (column 1) to plot}
  \item{ymin}{Minimum y-value (column 2) to plot}
  \item{ymax}{Maximum y-value (column 2) to plot}  
  \item{plotype}{Type of plot to generate: 1= points and lines, 2 = points, 3 = lines}
  \item{genplot}{Generate summary plots? (T or F)}
  \item{verbose}{Verbose output? (T or F)}
}

\seealso{ \code{\link{idPts}}, \code{\link{iso}}, \code{\link{trim}} and \code{\link{trimAT}} }