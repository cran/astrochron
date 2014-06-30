\name{autoPlot}
\alias{autoPlot}
\title{Automatically plot stratigraphic series, with smoothing if desired}
\description{
Automatically plot and smooth specified stratigraphic data, versus location. Data are smoothed with a Gaussian kernel.
}
\usage{
autoPlot(dat,cols=NULL,nrows=NULL,smooth=0,xgrid=1,output=F,verbose=T)
}
\arguments{
  \item{dat}{Your data frame; first column should be location identifier (e.g., depth).}
  \item{cols}{A vector that identifies the columns to extract (first column automatically extracted).}
  \item{nrows}{Number of rows in figure.}
  \item{smooth}{Width (temporal or spatial dimension) for smoothing with a Gaussian kernel (0 = no smoothing); the Gaussian kernel is scaled so that its quartiles (viewed as probability densities, that is, containing 50 percent of the area) are at +/- 25 percent of this value.}
  \item{xgrid}{For kernal smoothing: (1) evaluate on ORIGINAL sample grid, or (2) evaluate on EVENLY SPACED grid covering range.}
  \item{output}{Output data frame of smoothed values? (T or F)}
  \item{verbose}{verbose output (T or F)}
}