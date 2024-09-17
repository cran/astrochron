\name{mtm}
\alias{mtm}
\title{Multitaper method spectral analysis}
\description{
Multitaper method (MTM) spectral analysis (Thomson, 1982)
}
\usage{
mtm(dat,tbw=3,ntap=NULL,padfac=5,demean=T,detrend=F,siglevel=0.9,ar1=T,output=0,
     CLpwr=T,xmin,xmax,pl=1,sigID=T,genplot=T,verbose=T)
}
\arguments{
  \item{dat}{Stratigraphic series for MTM spectral analysis. First column should be location (e.g., depth), second column should be data value.}
  \item{tbw}{MTM time-bandwidth product.}
  \item{ntap}{Number of DPSS tapers to use. By default, this is set to (2*tbw)-1.}
  \item{padfac}{Pad with zeros to (padfac*npts) points, where npts is the original number of data points.}
  \item{demean}{Remove mean from data series? (T or F)}
  \item{detrend}{Remove linear trend from data series? (T or F)}
  \item{siglevel}{Significance level for peak identification. (0-1)}
  \item{ar1}{Estimate conventional AR(1) noise spectrum and confidence levels? (T or F)}
  \item{CLpwr}{Plot AR(1) noise confidence levels on power spectrum? (T or F)}
  \item{output}{What should be returned as a data frame? (0=nothing; 1= power spectrum + harmonic CL + AR1 CL + AR1 fit + 90\%-99\% AR1 power levels (ar1 must be set to TRUE to output AR model results); 2=significant peak frequencies; 3=significant peak frequencies + harmonic CL; 4=internal variables from spec.mtm). Option 4 is intended for expert users, and should generally be avoided.}
  \item{xmin}{Smallest frequency for plotting.}
  \item{xmax}{Largest frequency for plotting.}
  \item{pl}{Power spectrum plotting: (1) linear frequency-log spectral power, (2) linear frequency-linear spectral power (3) log frequency-log spectral power, (4) log frequency-linear spectral power}
  \item{sigID}{Identify significant frequencies on power and probabilty plots? (T or F)}
  \item{genplot}{Generate summary plots? (T or F)}
  \item{verbose}{Verbose output? (T or F)}
}

\details{
The power spectrum normalization approach applied here divides the Fourier coefficients by the 
number of points (npts) in the stratigraphic series, which is equivalent to dividing the power 
by (npts*npts). The (npts*npts) normalization has the convenient property whereby
-- for an unpadded series -- the sum of the power in the positive frequencies is equivalent 
to half of variance; the other half of the variance is in the negative frequencies.
 
If ar1=T, candidiate astronomical cycles are identified via isolation of those frequencies 
that achieve the required (e.g., 90 percent) "red noise" confidence level and MTM harmonic F-test 
confidence level. Allowance is made for the smoothing inherent in the MTM power spectral estimate 
as compared to the MTM harmonic spectrum. That is, an F-test peak is reported if it achieves the 
required MTM harmonic confidence level, while also achieving the required red noise confidence level 
within +/- half the power spectrum bandwidth resolution. One additional criterion is included to 
further reduce the false positive rate, a requirement that significant F-tests must occur on a local 
power spectrum high, which is parameterized as occurring above the local red noise background estimate.
See Meyers (2012) for futher information.

Note that the 'spec.mtm' function in package 'multitaper' (Rahim et al., 2014) is used for MTM spectrum estimation.

}

\references{
S.R. Meyers, 2012,
\emph{Seeing Red in Cyclic Stratigraphy: Spectral Noise Estimation for Astrochronology}: Paleoceanography, 27, PA3228, doi:10.1029/2012PA002307.

Rahim, K.J. and Burr W.S. and Thomson, D.J., 2014, \emph{Applications of Multitaper Spectral Analysis to Nonstationary Data}. PhD thesis, Queen's University. R package version 1.0-17, https://CRAN.R-project.org/package=multitaper.

Thomson, D. J., 1982, \emph{Spectrum estimation and harmonic analysis}, Proc. IEEE, 70, 1055-1096, doi:10.1109/PROC.1982.12433.
}

\seealso{ \code{\link{eha}}, \code{\link{lowspec}}, \code{\link{mtmAR}}, \code{\link{mtmML96}}, \code{\link{periodogram}}}

\examples{
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