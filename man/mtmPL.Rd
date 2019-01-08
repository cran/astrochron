\name{mtmPL}
\alias{mtmPL}
\title{Multitaper Method Spectral Analysis with Power Law (1/f) fit}
\description{
Multitaper Method (MTM) Spectral Analysis with Power Law (1/f) fit
}
\usage{
mtmPL(dat,tbw=3,ntap=NULL,padfac=5,demean=T,detrend=F,siglevel=0.9,flow=NULL,fhigh=NULL,
    output=0,CLpwr=T,xmin=0,xmax=Nyq,pl=1,sigID=F,genplot=T,verbose=T)
}
\arguments{
  \item{dat}{Stratigraphic series for MTM spectral analysis. First column should be location (e.g., depth), second column should be data value.}
  \item{tbw}{MTM time-bandwidth product.}
  \item{ntap}{Number of DPSS tapers to use. By default, this is set to (2*tbw)-1.}
  \item{padfac}{Pad with zeros to (padfac*npts) points, where npts is the original number of data points.}
  \item{demean}{Remove mean from data series? (T or F)}
  \item{detrend}{Remove linear trend from data series? (T or F)}
  \item{siglevel}{Significance level for peak identification.}
  \item{flow}{Lowest frequency to include in 1/f fit}
  \item{fhigh}{Highest frequency to include in 1/f fit}
  \item{output}{What should be returned as a data frame? (0=nothing; 1=spectrum + CLs + power law fit; 2=sig peak freqs; 3=sig peak freqs + prob; 4=all)}
  \item{CLpwr}{Plot power law noise confidence levels on power spectrum (in addition to the power law fit)? (T or F)}
  \item{xmin}{Smallest frequency for plotting.}
  \item{xmax}{Largest frequency for plotting.}
  \item{pl}{Power spectrum plotting: (1) linear frequency-log spectral power, (2) linear frequency-linear spectral power (3) log frequency-log spectral power, (4) log frequency-linear spectral power}
  \item{sigID}{Identify signficant frequencies on power and probabilty plots? (T or F)}
  \item{genplot}{Generate summary plots? (T or F)}
  \item{verbose}{Verbose output? (T or F)}
}

\details{
Candidiate astronomical cycles are identified via isolation of those frequencies 
that achieve the required (e.g., 90 percent) power law confidence level and MTM harmonic F-test 
confidence level. Allowance is made for the smoothing inherent in the MTM power spectral estimate 
as compared to the MTM harmonic spectrum. That is, an F-test peak is reported if it achieves the 
required MTM harmonic confidence level, while also achieving the required power law confidence level 
within +/- half the power spectrum bandwidth resolution. One additional criterion is included to 
further reduce the false positive rate, a requirement that significant F-tests must occur on a local 
power spectrum high, which is parameterized as occurring above the local red noise background estimate.
See Meyers (2012) for futher information.
}

\references{
S.R. Meyers, 2012,
\emph{Seeing Red in Cyclic Stratigraphy: Spectral Noise Estimation for Astrochronology}: Paleoceanography, 27, PA3228, doi:10.1029/2012PA002307.

Rahim, K.J. and Burr W.S. and Thomson, D.J., 2014, 
\emph{Appendix A: Multitaper R package in "Applications of Multitaper Spectral Analysis to Nonstationary Data"}, PhD diss., Queen's Univieristy, pp 149-183. http://hdl.handle.net/1974/12584

Thomson, D. J., 1982, 
\emph{Spectrum estimation and harmonic analysis}, Proc. IEEE, 70, 1055-1096, doi:10.1109/PROC.1982.12433.
}

\seealso{ \code{\link{eha}}, \code{\link{lowspec}}, \code{\link{mtm}}, \code{\link{mtmAR}}, \code{\link{mtmML96}}, \code{\link{periodogram}}, and \code{\link{spec.mtm}}}