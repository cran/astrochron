\name{mtmAR}
\alias{mtmAR}
\title{Intermediate spectrum test of Thomson et al. (2001)}
\description{
Perform the 'intermediate spectrum test' of Thomson et al. (2001).

Paraphrased from Thomson et al. (2001): 
Form an intermediate spectrum by dividing MTM by AR estimate. 
Choose an order P for a predictor. A variety of formal methods are available in the literature, 
but practically, one keeps increasing P (the order) until the range of the intermediate spectrum 
Si(f) (equation (C4) of Thomson et al., 2001) stops decreasing rapidly as a function of P. If the 
intermediate spectrum is not roughly white, as judged by the minima, the value of P should be 
increased. 
}
\usage{
mtmAR(dat,tbw=3,ntap=NULL,order=1,method="mle",CItype=1,padfac=5,demean=T,detrend=F,
      output=1,xmin=0,xmax=Nyq,pl=1,genplot=T,verbose=T)
}
\arguments{
  \item{dat}{Stratigraphic series for analysis. First column should be location (e.g., depth), second column should be data value.}
  \item{tbw}{MTM time-bandwidth product.}
  \item{ntap}{Number of DPSS tapers to use. By default, this is set to (2*tbw)-1.}
  \item{order}{Order of the AR spectrum.}
  \item{method}{AR method ("yule-walker", "burg", "ols", "mle", "yw")}
  \item{CItype}{Illustrate (1) one-sided or (2) two-sided confidence intervals on plots}
  \item{padfac}{Pad with zeros to (padfac*npts) points, where npts is the original number of data points.}
  \item{demean}{Remove mean from data series? (T or F)}
  \item{detrend}{Remove linear trend from data series? (T or F)}
  \item{output}{Output (1) intermediate spectrum and confidence levels, (2) intermediate spectrum, (3) confidence levels}
  \item{xmin}{Smallest frequency for plotting.}
  \item{xmax}{Largest frequency for plotting.}
  \item{pl}{Plot logarithm of spectral power (1) or linear spectral power (2)?}
  \item{genplot}{Generate summary plots? (T or F)}
  \item{verbose}{Verbose output? (T or F)}
}

\details{
The power spectrum normalization approach applied here divides the Fourier coefficients by the 
number of points (npts) in the stratigraphic series, which is equivalent to dividing the power 
by (npts*npts). The (npts*npts) normalization has the convenient property whereby
-- for an unpadded series -- the sum of the power in the positive frequencies is equivalent 
to half of variance; the other half of the variance is in the negative frequencies.

Note that the 'spec.mtm' function in package 'multitaper' (Rahim et al., 2014) is used for MTM spectrum estimation.

}

\references{
Rahim, K.J. and Burr W.S. and Thomson, D.J., 2014, \emph{Applications of Multitaper Spectral Analysis to Nonstationary Data}. PhD thesis, Queen's University. R package version 1.0-17, https://CRAN.R-project.org/package=multitaper.

Thomson, D. J., L. J. Lanzerotti, and C. G. Maclennan, 2001, \emph{The interplanetary magnetic field: Statistical properties and discrete modes}, J. Geophys.Res., 106, 15,941-15,962, doi:10.1029/2000JA000113.
}

\seealso{ \code{\link{eha}}, \code{\link{lowspec}}, \code{\link{mtm}}, \code{\link{mtmML96}}, and \code{\link{periodogram}}}

\examples{
# generate example series with periods of 400 ka, 100 ka, 40 ka and 20 ka
ex = cycles(freqs=c(1/400,1/100,1/40,1/20),start=1,end=1000,dt=5)

# add AR1 noise
noise = ar1(npts=200,dt=5,sd=.5)
ex[2] = ex[2] + noise[2]

# MTM spectral analysis, with conventional AR1 noise test
pl(1,title="mtmAR")
mtmAR(ex)
}