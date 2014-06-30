\name{zoomIn}
\alias{zoomIn}
\title{Dynamically explore cross-plot, zoom-in into specified region}
\description{
Dynamically explore cross-plot, zoom-in into specfied region. Accepts one dataframe/matrix 
with two columns, or two dataframes/vectors with one column.
}
\usage{
zoomIn(dat1,dat2=NULL,ptsize=1,xmin=NULL,xmax=NULL,ymin=NULL,ymax=NULL,plotype=1,
       verbose=T)
}
\arguments{
  \item{dat1}{Data frame with one or two columns. If one column, dat2 must also be specified.}
  \item{dat2}{Data frame with one column.}
  \item{ptsize}{Size of plotted points.}
  \item{xmin}{Minimum x-value (column 1) to plot}
  \item{xmax}{Maximum x-value (column 1) to plot}  
  \item{ymin}{Minimum y-value (column 2) to plot}
  \item{ymax}{Maximum y-value (column 2) to plot}  
  \item{plotype}{Type of plot to generate: 1= points and lines, 2 = points, 3 = lines}
  \item{verbose}{Verbose output? (T or F)}
}