\name{delPts}
\alias{delPts}
\title{Interactively Delete Points in Plot}
\description{
Interactively delete points in x,y plot.
}
\usage{
delPts(dat,ptsize=1,xmin=NULL,xmax=NULL,ymin=NULL,ymax=NULL,plotype=1)
}
\arguments{
  \item{dat}{Data frame with two columns}
  \item{ptsize}{Size of plotted points.}
  \item{xmin}{Minimum x-value (column 1) to plot}
  \item{xmax}{Maximum x-value (column 1) to plot}
  \item{ymin}{Minimum y-value (column 2) to plot}
  \item{ymax}{Maximum y-value (column 2) to plot}  
  \item{plotype}{Type of plot to generate: 1= points and lines, 2 = points, 3 = lines}
}

\seealso{ \code{\link{idPts}}, \code{\link{iso}}, \code{\link{trim}} and \code{\link{trimAT}} }