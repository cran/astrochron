\name{slideCor}
\alias{slideCor}
\title{Identify optimal spatial/temporal shift to maximize correlation between two stratigraphic/time series.}
\description{
Identify optimal spatial/temporal shift to maximize correlation between two stratigraphic/time series.
}
\usage{
slideCor(dat1,dat2,rev=F,cormethod=1,minpts=NULL,detrend=F,rmin=NULL,
         output=T,genplot=T,verbose=T)
}
\arguments{
  \item{dat1}{Stratigraphic series 1. First column should be location (e.g., depth), second column should be data value.}
  \item{dat2}{Stratigraphic series 2. First column should be location (e.g., depth), second column should be data value.}
  \item{rev}{Reverse polarity of stratigraphic series 2 (multiply proxy data value by -1)? (T or F)}
  \item{cormethod}{Method used for calculation of correlation coefficient (1=Pearson, 2=Spearman rank, 3=Kendall)}
  \item{minpts}{Minimum number of data points for calculation of correlation coefficient.}
  \item{detrend}{Remove linear trend from each window? (T or F)}
  \item{rmin}{Minimum r and r2 value shown on plots. By default all r and r2 values will be displayed.}
  \item{output}{Output correlation coefficient results as a dataframe? (T or F)}
  \item{genplot}{Generate summary plots? (T or F)}
  \item{verbose}{Verbose output? (T or F)}
}

\details{
slideCor is a general purpose tool to identify the optimal spatial/temporal correlation 
between two data sets. A few example applications include:
(1) stratigraphic correlation of data series from two locations (as in Preto et al., 2004), 
(2) identification of the optimal spatial/temporal lag between two variables from the 
same site, and (3) identification of the optimal fit between a floating astrochronology 
and astronomical target (e.g, Mitchell et al., 2008). 

Both series must be evenly sampled, but are not required to have the same sampling interval.
If stratigraphic series of different duration/length are being compared, the shift (in 
spatial or temporal units) should be interpreted as the location within the longer
stratigraphic series where the shorter stratigraphic series begins.
If both stratigraphic series are of the same duration/length, then the shift is the
location within dat1 where dat2 begins.

In some cases, it may be desirable to smooth or bandpass the data series before 
implementing slideCor (e.g., functions noLow, noKernel, bandpass, taner, etc.).
}

\references{
Preto, N., Hinnov, L.A., De Zanche, V., Mietto, P., and Hardie, L.A., 2004, \emph{The Milankovitch interpretation of the Latemar Platform cycles (Dolomites, Italy): Implications for geochronology, biostratigraphy, and Middle Triassic carbonate accumulation}, SEPM Special Publication 81.

Mitchell, R.N., Bice, D.M., Montanari, A., Cleaveland, L.C., Christianson, K.T., Coccioni, R., and Hinnov, L.A., 2008, \emph{Oceanic anoxic cycles? Orbital prelude to the Bonarelli Level (OAE 2)}, Earth Planet. Sci. Lett. 26, 1-16.
}


\seealso{ \code{\link{surrogateCor}} }

\examples{
# Example 1: generate AR1 noise
ex1 <- ar1(npts=1000,dt=1)
# isolate a section
ex2 <- iso(ex1,xmin=200,500)
ex2[1] <- ex2[1]-200

res=slideCor(ex1,ex2)

# Example 2: an astronomical signal
ex1=etp(tmin=0,tmax=1000)
# isolate a 200 ka section
ex2=iso(ex1,xmin=400,xmax=600)
# convert to a floating timescale (elapsed time)
ex2[1] <- ex2[1]-400

res=slideCor(ex1,ex2)
# now anchor the floating time scale
anchor <- ex1[res[which.max(res[,2]),1],1]
ex2.anchor <- anchorTime(dat=ex2, time=0, age=anchor, timeDir=2)
}