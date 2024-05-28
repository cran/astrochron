\name{stratPCA}
\alias{stratPCA}
\title{Principal Component Analysis on Stratigraphic Series}
\description{
Perform principal component analysis on stratigraphic series.  
}
\usage{
stratPCA(dat,id=TRUE,rot=0,nPC=NULL,output=0,symSize=2,genplot=1)
}
\arguments{
  \item{dat}{Your data frame; should be in stratigraphic order. First column may be location identifier (e.g., depth). Can contain any number of variables.}
  \item{id}{Is the first column of dat an ID column (depth/height/time)? (T or F)}
  \item{rot}{Rotation method. 0=none, 1=varimax (orthogonal), 2=promax (oblique)}
  \item{nPC}{Number of principal components to extract. Default=all}
  \item{output}{What would you like to output? 0=nothing, 1=PC loadings, 2=scores}
  \item{symSize}{Size of symbols used in color plots? Default=2}
  \item{genplot}{Generate summary plots? (0) no plots, (1) standard plots, (2) additional plots}
}
\examples{
\dontrun{
# create a test data set, composed of 7 variables:
#  antiphased 20 kyr cycles with noise, 
#  antiphased 40 kyr cycles with noise,
#  antiphased 110 kyr cycles with noise, 
#  and a variable that is entirely noise
noise=0.01
a=cycles(1/20,noisevar=noise,genplot=FALSE)
b=cycles(1/20,phase=pi,noisevar=noise,genplot=FALSE)
c=cycles(1/40,noisevar=noise,genplot=FALSE)
d=cycles(1/40,phase=pi,noisevar=noise,genplot=FALSE)
e=cycles(1/110,noisevar=noise,genplot=FALSE)
f=cycles(1/110,phase=pi, noisevar=noise,genplot=FALSE)
g=ar1(npts=500,genplot=FALSE)

ex=data.frame(cbind(a[,1],a[,2],b[,2],c[,2],d[,2],e[,2],f[,2],g[,2]))

stratPCA(ex)

stratPCA(ex,rot=1,nPC=4)

stratPCA(ex,rot=2,nPC=4)
}
}