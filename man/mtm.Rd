\name{mtm}
\alias{mtm}
\title{Multitaper Method Spectral Analysis}
\description{
Multitaper Method (MTM) Spectral Analysis
}
\usage{
mtm(dat,tbw=3,ntap=NULL,padfac=5,demean=T,detrend=F,siglevel=0.9,ar1=T,output=0,
     CLpwr=T,xmin,xmax,pl=1,sigID=F,genplot=T,verbose=T)
}
\arguments{
  \item{dat}{Stratigraphic series for MTM spectral analysis. First column should be location (e.g., depth), second column should be data value.}
  \item{tbw}{MTM time-bandwidth product.}
  \item{ntap}{Number of DPSS tapers to use. By default, this is set to (2*tbw)-1.}
  \item{padfac}{Pad with zeros to (padfac*npts) points, where npts is the original number of data points.}
  \item{demean}{Remove mean from data series? (T or F)}
  \item{detrend}{Remove linear trend from data series? (T or F)}
  \item{siglevel}{Significance level for peak identification.}
  \item{ar1}{Estimate conventional AR(1) noise spectrum and confidence levels? (T or F)}
  \item{CLpwr}{Plot AR(1) noise confidence levels on power spectrum? (T or F)}
  \item{output}{What should be returned as a data frame? (0=nothing; 1=spectrum + CLs + AR1 fit; 2=sig peak freqs; 3=sig peak freqs + prob; 4=all)}
  \item{xmin}{Smallest frequency for plotting.}
  \item{xmax}{Largest frequency for plotting.}
  \item{pl}{Plot logarithm of spectral power (1) or linear spectral power (2)?}
  \item{sigID}{Identify signficant frequencies on power and probabilty plots? (T or F)}
  \item{genplot}{Generate summary plots? (T or F)}
  \item{verbose}{Verbose output? (T or F)}
}
\references{
Rahim, K.J. and Burr W.S. and Thomson, D.J., 2014, \emph{Appendix A: Multitaper R package in "Applications of Multitaper Spectral Analysis to Nonstationary Data"}, PhD diss., Queen's Univieristy, pp 149-183. http://hdl.handle.net/1974/12584

Thomson, D. J., 1982, \emph{Spectrum estimation and harmonic analysis}, Proc. IEEE, 70, 1055-1096, doi:10.1109/PROC.1982.12433.
}

\seealso{ \code{\link{spec.mtm}}, \code{\link{lowspec}}, and \code{\link{periodogram}} }

\examples{
\donttest{
# generate example series with periods of 400 ka, 100 ka, 40 ka and 20 ka
ex = cycles(freqs=c(1/400,1/100,1/40,1/20),start=1,end=1000,dt=5)

# add AR1 noise
noise = ar1(npts=200,dt=5,sd=.5)
ex[2] = ex[2] + noise[2]

# MTM spectral analysis, with conventional AR1 noise test
pl(1,title="mtm")
mtm(ex,ar1=TRUE)

# compare to ML96 analysis
pl(1, title="mtmML96")
mtmML96(ex)

# compare to analysis with LOWSPEC
pl(1, title="lowspec")
lowspec(ex)

# compare to amplitudes from eha
pl(1,title="eha")
eha(ex,tbw=3,win=1000,pad=1000)

}

}