\name{freq2sedrate}
\alias{freq2sedrate}
\title{Convert record of local spatial frequency (from EHA) to sedimentation rate curve}
\description{
Convert record of local spatial frequency (from EHA) to sedimentation rate curve
}
\usage{
freq2sedrate(freqs,period=NULL,ydir=1,genplot=T,verbose=T)
}
\arguments{
  \item{freqs}{Data frame containing depth/height in first column (meters) and spatial frequencies in second column (cycles/m)}
  \item{period}{Temporal period of spatial frequency (ka)}
  \item{ydir}{Direction for y-axis in plots (depth,height). -1 = values increase downwards (slower), 1 = values increase upwards }
  \item{genplot}{Generate summary plots? (T or F)}
  \item{verbose}{Verbose output? (T or F)}
}