\name{autoPlot}
\alias{autoPlot}
\title{Automatically plot multiple stratigraphic series, with smoothing if desired}
\description{
Automatically plot and smooth specified stratigraphic data, versus location. Data are smoothed with a Gaussian kernel if desired.
}
\usage{
autoPlot(dat,cols=NULL,vertical=T,ydir=NULL,nrows=NULL,plotype=1,smooth=0,xgrid=1,
        output=F,genplot=T,verbose=T)
}
\arguments{
  \item{dat}{Your data frame; first column should be location identifier (e.g., depth).}
  \item{cols}{A vector that identifies the columns to extract (first column automatically extracted).}
  \item{vertical}{Generate vertical stratigraphic plots? (T or F) If F, will generate horizontal plots.}
  \item{ydir}{Direction for stratigraphic axis in plots (depth,height,time). If vertical=T, then -1 results in values increasing downwards, while 1 results in values increasing upwards. If vertical=F, then -1 results in values increasing toward the left, while 1 results in values increasing toward the right.}
  \item{nrows}{Number of rows in figure (if vertical = T; otherwise this will be the number of columns).}
  \item{plotype}{Type of plot to generate: 1= points and lines, 2 = points, 3 = lines}
  \item{smooth}{Width (temporal or spatial dimension) for smoothing with a Gaussian kernel (0 = no smoothing); the Gaussian kernel is scaled so that its quartiles (viewed as probability densities, that is, containing 50 percent of the area) are at +/- 25 percent of this value.}
  \item{xgrid}{For kernal smoothing: (1) evaluate on ORIGINAL sample grid, or (2) evaluate on EVENLY SPACED grid covering range.}
  \item{output}{Output data frame of smoothed values? (T or F)}
  \item{genplot}{Generate summary plots (T or F)}
  \item{verbose}{Verbose output (T or F)}
}