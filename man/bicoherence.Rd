\name{bicoherence}
\alias{bicoherence}
\title{bicoherence: Calculate bispectrum and bicoherence using WOSA method as detailed in Choudhury et al. (2008).}
\description{
bicoherence: Calculate bispectrum and bicoherence using WOSA method as detailed in Choudhury et al. (2008). Kim and Powers (1979) bicoherence normalization is used, with a conditioning factor to reduce the potential for spurious bicoherence peaks.
}
\usage{
bicoherence(dat,overlap=50,segments=8,CF=0,CL=95,padfac=2,demean=T,detrend=T,
            taper=T,maxF=Nyq,output=0,genplot=T,color=1,id=NULL,logpwr=F,logbis=F,
            check=T,verbose=T)
}
\arguments{
  \item{dat}{Stratigraphic series to analyze. First column should be location (e.g., depth, time), second column should be data value.}
  \item{overlap}{Percent overlap for WOSA segments (use 0-50 percent).}
  \item{segments}{Number of segments for WOSA.}
  \item{CF}{Conditioning factor for bicoherence estimation (see pg. 70 of Choudhury et al., 2008). When CF=0, this will be automatically estimated as the 75th percentile of the magnitude-squared bicoherence denominator}
  \item{CL}{Confidence level to identify with a contour on the plots (0-100).}
  \item{padfac}{Pad with zeros to (padfac*npts) points, where npts is the original number of data points per segment.}
  \item{demean}{Remove mean from data series? (T or F)}
  \item{detrend}{Remove linear trend from data series? (T or F)}
  \item{taper}{Apply Hanning taper? (T or F)}
  \item{maxF}{Maximum frequency for analysis/plotting.}
  \item{output}{Return output as new data frame? 0 = none, 1 = magnitude-squared bispectrum, 2 = magnitude-squared bicoherence, 3 = magnitude-squared bicoherence confidence levels, 4 = everything.}
  \item{genplot}{Generate summary plots? (T or F)}
  \item{color}{Use (1) grayscale or (2) viridis color scale for 3D plots?}
  \item{id}{A vector listing frequencies to identify on the plots as diagonal lines.}
  \item{logpwr}{Use a log scale for power spectrum? (T or F)}
  \item{logbis}{Use a log scale for bispectrum? (T or F)}
  \item{check}{Conduct compliance checks before processing? (T or F) In general this should be activated; the option is included for Monte Carlo simulation.}
  \item{verbose}{Verbose output? (T or F)}
}

\details{
This function accompanies the publication Sullivan et al. (2023, PNAS): "Bicoherence is a higher-order statistic, which 
quantifies coupling of individual frequencies (as combination and difference tones). For three frequencies, f1, f2 and f3, 
high bicoherence indicates that f1 + f2 = f3, with the additional property that their phases are also coupled in the data 
series. We employ the bicoherence approach of Choudhury, Shah and Thornhill (2008), which uses Welch Overlapping Spectral 
Analysis (WOSA), as this method includes improvements to the standard approach that reduce the potential for spurious 
bicoherence peaks. Note that reliable bicoherence estimates require a substantial increase in bandwidth resolution ("smoothing")".
}

\references{
Choudhury, A.A.S., Shah, S.L., and Thornhill, N.F. (2008),
\emph{Diagnosis of Process Nonlinearities and Valve Stiction: Data Driven Approaches}: Springer, 284 pp.

Kim, Y.C., and Powers, E.J. (1979), 
\emph{Digital Bispectral Analysis and Its Applications to Nonlinear Wave Interactions}: IEEE Transactions on Plasma Science, v. PS-7, 120-131.

Sullivan, N.B., Meyers, S.R., Levy, R.H., McKay, R.M., Golledge, N.R., Cortese, G. (2023),
\emph{Millennial-scale variability of the Antarctic Ice Sheet during the Early Miocene}: Proceedings of the National Academy of Sciences, https://doi.org/10.1073/pnas.2304152120.

Welch, P.D. (1967),
\emph{The use of Fast Fourier Tranform for the estimation of power spectra: A method based on time averaging over short, modified periodograms}: IEEE Transactions on Audio and Electroacoustics, AU-5 (2): 70-73.
}

\examples{
 \dontrun{
# Generate an example quadratic phase-coupled signal as in Choudhury et al. (2008, pg. 79)
    n = rnorm(500)
    t=1:500
    signal = sin(2*pi*0.12*t + pi/3) + sin(2*pi*0.18*t + pi/12) + sin(2*pi*0.3*t + 5*pi/12) + n
    ex=data.frame(cbind(t,signal))
    bicoherence(ex)
 }
}