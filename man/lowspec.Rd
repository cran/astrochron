\name{lowspec}
\alias{lowspec}
\title{Robust Locally-Weighted Regression Spectral Background Estimation}
\description{
LOWSPEC: Robust Locally-Weighted Regression Spectral Background Estimation
}
\usage{
lowspec(dat,decimate=NULL,tbw=3,padfac=5,detrend=F,siglevel=0.9,xmin,xmax,
        setrho,lowspan,b_tun,output=0,sigID=T,pl=1,genplot=T,verbose=T)
}
\arguments{
  \item{dat}{Stratigraphic series for LOWSPEC. First column should be location (e.g., depth), second column should be data value.}
  \item{decimate}{Decimate statigraphic series to have this sampling interval (via piecewise linear interpolation). By default, no decimation is performed.}
  \item{tbw}{MTM time-bandwidth product (2 or 3 permitted)}
  \item{padfac}{Pad with zeros to (padfac*npts) points, where npts is the original number of data points.}
  \item{detrend}{Remove linear trend from data series? This detrending is performed following AR1 prewhitening. (T or F)}
  \item{siglevel}{Significance level for peak identification.}
  \item{xmin}{Smallest frequency for plotting.}
  \item{xmax}{Largest frequency for plotting.}
  \item{setrho}{Define AR1 coefficient if desired (otherwise calculated).}
  \item{lowspan}{Span for LOWESS smoothing of prewhitened signal, usually fixed to 1. If using value <1, method is overly conservative with reduced false positive rate.}
  \item{b_tun}{Robustness weight parameter for LOWSPEC.}
  \item{output}{What should be returned as a data frame? (0=nothing; 1=spectrum + CLs; 2=sig peaks; 3=PDF image)}
  \item{sigID}{Identify signficant frequencies on power and probabilty plots? (T or F)}
  \item{pl}{Plot logarithm of spectral power (1) or linear spectral power (2)?}
  \item{genplot}{Generate summary plots? (T or F)}
  \item{verbose}{Verbose output? (T or F)}
}

\details{ 
LOWSPEC is a 'robust' method for spectral background estimation, designed for the 
identification of potential astronomical signals that are imbedded in red noise (Meyers, 2012). 
The complete algoritm implemented here is as follows: (1) initial pre-whitening with AR1 
filter (default) or other filter as appropriate (see function prewhiteAR), (2) power spectral 
estimation via the multitaper method (Thomson, 1982), (3) robust locally weighted estimation of 
the spectral background using the LOWESS-based (Cleveland, 1979) procedure of Ruckstuhl et al. (2001), 
(4) assignment of confidence levels using a Chi-square distribution. 

Candidiate astronomical cycles are subsequently idenitified via isolation of those frequencies 
that achieve the required (e.g., 90 percent) LOWSPEC confidence level and MTM harmonic F test 
confidence level. Allowance is made for the smoothing inherent in the MTM power spectral estimate 
as compared to the MTM harmonic spectrum. That is, an F test peak is reported if it achieves the 
required MTM harmonic confidence level, while also achieving the required LOWSPEC confidence level 
within +/- half the power spectrum bandwidth resolution. One additional criterion is included to 
further reduce the false positive rate, a requirement that significant F tests must occur on a local 
power spectrum high, which is parameterized as occurring above the local LOWSPEC background estimate.
See Meyers (2012) for futher information on the algorithm.

In this implementation, the 'robustness criterion' ('b' in EQ. 6 of Ruckstuhl et al., 2001) has been 
optimized for 2 and 3 pi DPSS, using a 'span' of 1. By default the robustness criterion will be 
estimated. Both 'b' and the 'span' can be expliclty set using parameters 'b_tun' and 'lowspan'.
Note that it is permissible to decrease 'lowspan' from its default value, but this will result in 
an overly conservative false positive rate. However, it may be necessary to reduce 'lowspan' to provide 
an approporiate background fit for some stratigraphic data. Another option is to decimate the data
series prior to spectral estimation.
}

\value{
If option 1 is selected, a data frame containing the following is returned: Frequency, Prewhitened power, LOWSPEC background, LOWSPEC CL, F-test CL.

If option 2 is selected, the 'significant' frequencies are returned (as described above).

If option 3 is selected, the graphics are output to a PDF file.
}

\references{
W.S. Cleveland, 1979, 
\emph{Locally weighted regression and smoothing scatterplots}: Journal of the American Statistical Association, v. 74, p. 829-836.

S.R. Meyers, 2012,
\emph{Seeing Red in Cyclic Stratigraphy: Spectral Noise Estimation for Astrochronology}: Paleoceanography, 27, PA3228, doi:10.1029/2012PA002307.

A.F. Ruckstuhl, M.P Jacobson, R.W. Field, and J.A. Dodd, 2001,
\emph{Baseline subtraction using robust local regression estimation}: Journal of Quantitative Spectroscopy & Radiative Transfer, v. 68, p. 179-193.

D.J. Thomson, 1982, 
\emph{Spectrum estimation and harmonic analysis}: IEEE Proceedings, v. 70, p. 1055-1096.
}

\seealso{ \code{\link{spec.mtm}}, \code{\link{baseline}}, \code{\link{mtm}}, and \code{\link{periodogram}} }


\examples{
# generate example series with periods of 400 ka, 100 ka, 40 ka and 20 ka
ex = cycles(freqs=c(1/400,1/100,1/40,1/20),start=1,end=1000,dt=5)

# add AR1 noise
noise = ar1(npts=200,dt=5,sd=.5)
ex[2] = ex[2] + noise[2]

# LOWSPEC analysis
#lowspec(ex)
}