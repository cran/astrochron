\name{mtmML96}
\alias{mtmML96}
\title{Mann and Lees (1996) robust red noise MTM analysis}
\description{
Mann and Lees (1996) robust red noise MTM analysis. This function implements several improvements to the algorithm used in SSA-MTM toolkit, including faster AR1 model optimization, and more appropriate 'edge-effect' treatment.
}
\usage{
mtmML96(dat,tbw=3,ntap=NULL,padfac=5,demean=T,detrend=F,medsmooth=0.2,
         opt=1,linLog=2,siglevel=0.9,output=0,CLpwr=T,xmin=0,xmax=Nyq,
         sigID=F,pl=1,genplot=T,verbose=T)
}
\arguments{
  \item{dat}{Stratigraphic series for MTM spectral analysis. First column should be location (e.g., depth), second column should be data value.}
  \item{tbw}{MTM time-bandwidth product.}
  \item{ntap}{Number of DPSS tapers to use. By default, this is set to (2*tbw)-1.}
  \item{padfac}{Pad with zeros to (padfac*npts) points, where npts is the original number of data points.}
  \item{demean}{Remove mean from data series? (T or F)}
  \item{detrend}{Remove linear trend from data series? (T or F)}
  \item{medsmooth}{ML96 median smoothing parameter (1 = use 100 percent of spectrum; 0.20 = use 20 percent)}
  \item{opt}{Optimization method for robust AR1 model estimation (1=Brent's method:fast, 2=Gauss-Newton:fast, 3=grid search:slow)}
  \item{linLog}{Optimize AR1 model fit using (1) linear power or (2) log(power)?}
  \item{siglevel}{Significance level for peak identification. (0-1)}
  \item{output}{What should be returned as a data frame? (0=nothing; 1=spectrum + CLs + AR1 fit + median smoothed spectrum; 2=significant peak frequencies; 3=significant peak frequencies + harmonic CL)}
  \item{CLpwr}{Plot ML96 AR(1) noise confidence levels on power spectrum? (T or F)}
  \item{xmin}{Smallest frequency for plotting.}
  \item{xmax}{Largest frequency for plotting.}
  \item{sigID}{Identify significant frequencies on power and probabilty plots? (T or F)}
  \item{pl}{Plot logarithm of spectral power (1) or linear spectral power (2)?}
  \item{genplot}{Generate summary plots? (T or F)}
  \item{verbose}{Verbose output? (T or F)}
}

\details{
This function conducts the Mann and Lees (1996; ML96) "robust red noise" analysis, with an improved median
smoothing approach. The original Mann and Lees (1996) approach applies a truncation of the median 
smoothing window to include fewer frequencies near the edges of the spectrum; while truncation is
required, its implementation in the original method often results in an "edge effect" that can produce 
excess false positive rates at low frequencies, commonly within the eccentricity-band (Meyers, 2012). 

To help address this issue, an alternative median smoothing approach is applied that implements Tukey's 
robust end-point rule and symmetrical medians (see the function runmed for details). Numerical experiments 
indicate that this approach produces an approximately uniform false positive rate across the spectrum.  
It should be noted that the false positive rates are still inflated with this method, but they are substantially
reduced compared to the original ML96 approach. For example, simulations using rho=0.9 (using identical 
parameters to those in Meyers, 2012) yield median false positive rates of 1.7\%, 7.3\% and 13.4\%, for the 99\%, 
95\% and 90\% confidence levels (respectively). This compares with 4.7\%, 11.4\% and 17.8\% using 
the original approach (see Table 2 of Meyers, 2012).

NOTE: If the (fast) Brent or Gauss-Newton methods fail, use the (slow) grid search approach.

This version of the ML96 algorithm was first implemented in Patterson et al. (2014).
}

\references{
Mann, M.E., and Lees, J.M., 1996, \emph{Robust estimation of background noise and signal 
detection in climatic time series}, Clim. Change, 33, 409-445.

Meyers, S.R., 2012, \emph{Seeing red in cyclic stratigraphy: Spectral noise estimation for
astrochronology}, Paleoceanography, 27, PA3228.

Patterson, M.O., McKay, R., Naish, T., Escutia, C., Jimenez-Espejo, F.J., 
Raymo, M.E., Meyers, S.R., Tauxe, L., Brinkhuis, H., and IODP Expedition 318 Scientists,
2014, \emph{Response of the East Antarctic Ice Sheet to orbital forcing during the 
Pliocene and Early Pleistocene}, Nature Geoscience, v. 7, p. 841-847.

Thomson, D. J., 1982, \emph{Spectrum estimation and harmonic analysis}, Proc. IEEE, 
70, 1055-1096, doi:10.1109/PROC.1982.12433.

http://www.meteo.psu.edu/holocene/public_html/Mann/tools/MTM-RED/

Tukey, J.W., 1977, \emph{Exploratory Data Analysis}, Addison.
}

\seealso{ \code{\link{eha}}, \code{\link{lowspec}}, \code{\link{mtm}}, \code{\link{mtmAR}}, \code{\link{periodogram}}, \code{\link{runmed}}, and \code{\link{spec.mtm}}}

\examples{
\donttest{
# generate example series with periods of 400 ka, 100 ka, 40 ka and 20 ka
ex = cycles(freqs=c(1/400,1/100,1/40,1/20),start=1,end=1000,dt=5)

# add AR1 noise
noise = ar1(npts=200,dt=5,sd=0.5)
ex[2] = ex[2] + noise[2]

# run ML96 analysis
pl(1, title="mtmML96")
mtmML96(ex)

# compare to analysis with conventional AR1 noise test
pl(1,title="mtm")
mtm(ex)

# compare to analysis with LOWSPEC
pl(1, title="lowspec")
lowspec(ex)

# compare to amplitudes from eha
pl(1,title="eha")
eha(ex,tbw=3,win=1000,pad=1000)

}
}
