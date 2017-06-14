\name{eAsmTrack}
\alias{eAsmTrack}
\title{Track ASM null hypothesis significance level minima in eASM results}
\description{
Track ASM null hypothesis significance level minima in eASM results.
}
\usage{
eAsmTrack(res,threshold=.5,ydir=-1,genplot=T,verbose=T)
}
\arguments{
  \item{res}{eAsm results. Must have the following format: column 
  1=sedimentation rate; remaining columns (2 to n)=Ho-SL; titles for columns 2 to 
  n must be the location (depth or height). Note that this format is ouput by function eAsm.}
  \item{threshold}{Threshold Ho-SL value for analysis and plotting.}
  \item{ydir}{Direction for y-axis in plots (depth or height). -1 = values increase downwards (slower plotting!), 1 = values increase upwards.}
  \item{genplot}{Generate summary plots? (T or F)}
  \item{verbose}{Verbose output? (T or F)}
}

\details{
Please see function eAsm for details.
}

\seealso{ \code{\link{asm}}, \code{\link{eAsm}}, and \code{\link{eha}} }
