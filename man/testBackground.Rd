\name{testBackground}
\alias{testBackground}
\title{Evaluate power spectrum false positive rates via Monte Carlo simulation}
\description{
This is a simulation tool to evaluate power spectrum false positive rates, the 
frequency distribution of the false positives, and the behavior of numerous
"multiple correction" procedures, for a range of background estimation 
approaches that are implemented in Astrochron. The tool can be used to conduct
surrogate analyses, alongside analysis of real data, to better understand
the suitability of particular background estimation approaches. The resulting 
simulations are similar to those presented in Figure 3 of Meyers (2012) and 
Crampton et al. (PNAS).
}
\usage{
testBackground(npts=1001,dt=5,noiseType="ar1",coeff=NULL,method="periodogramAR1",
               opt=NULL,demean=T,detrend=F,low=0,tbw=3,multi=F,iter=2000,output=F,
               verbose=T,genplot=F)
}
\arguments{
  \item{npts}{Number of points in simulated stratigraphic series (surrogates).}
  \item{dt}{Sampling interval for surrogates.}
  \item{noiseType}{Select "ar1" for AR1 noise surrogates, or "pwrLaw" for Power Law noise surrogates}
  \item{coeff}{AR1 coefficient (rho) or Power Law coefficient (beta) for surrogates.}
  \item{method}{Background estimation method: (1) "mtmAR1" (function mtm), (2) "mtmML96" (function mtmML96), (3) "lowspec" (function lowspec), (4) "mtmPL" (function mtmPL), (5) "periodogramPL" (function periodogram), (6) "periodogramAR1" (function periodogram)}
  \item{opt}{Method specific options. For mtmML96, this is medsmooth (see function mtmML96); for lowspec this is lowspan (see function lowspec); for periodogram this is percent cosine taper (see function cosTaper).}
  \item{demean}{Remove mean value from simulated surrogates? (T or F; this option does not apply to lowspec)}
  \item{detrend}{Remove linear trend from simulated surrogates? (T or F)}
  \item{low}{Remove long-term trend using a LOWESS smoother? Choose a value ranging from 0-1 (see function noLow). 0 = no long-term trend removal.}
  \item{tbw}{MTM time-bandwidth product. This option is ignored for methods 5 and 6.}
  \item{multi}{Evaluate a range of multiple-comparison tests too? (T or F)}
  \item{iter}{Number of iterations (surrogate series) for Monte Carlo simulation.}
  \item{output}{Output data frame? (T or F)} 
  \item{verbose}{Verbose output? (T or F)}
  \item{genplot}{Generate summary plots? (T or F)}
}
\details{
The Monte Carlo simulations can utilize AR1 or Power Law noise surrogates.
Background estimation approaches include conventional AR1, ML96, LOWSPEC 
and Power Law. The function also allows evaluation of common data detrending
approaches (linear trend removal, LOWESS trend removal).

Note that MTM-ML96 conducts the Mann and Lees (1996; ML96) "robust red noise" analysis, 
with an improved median smoothing approach. The original Mann and Lees (1996) approach 
applies a truncation of the median smoothing window to include fewer frequencies near 
the edges of the spectrum; while truncation is required, its implementation in the 
original method often results in an "edge effect" that can produce excess false positive 
rates at low frequencies, commonly within the eccentricity-band (Meyers, 2012). To help 
address this issue, an alternative median smoothing approach is applied that implements 
Tukey's robust end-point rule and symmetrical medians (see the function mtmML96 for more
details). This version of the ML96 algorithm was first implemented in Patterson et al. (2014).

See function multiTest for more information on the multiple comparison tests evaluated.
}

\references{
W.S. Cleveland, 1979, 
\emph{Locally weighted regression and smoothing scatterplots}: Journal of the American Statistical Association, v. 74, p. 829-836.

J.S. Campton, S.R. Meyers, R.A. Cooper, P.M Sadler, M. Foote, D. Harte, 2018, 
\emph{Pacing of Paleozoic macroevolutionary rates by Milankovitch grand cycles}: Proceedings of the National Academy of Sciences, doi:10.1073/pnas.1714342115.

M.E. Mann, and J.M. Lees, 1996,
 \emph{Robust estimation of background noise and signal detection in climatic time series}, Clim. Change, 33, 409-445.

S.R. Meyers, 2012,
\emph{Seeing Red in Cyclic Stratigraphy: Spectral Noise Estimation for Astrochronology}: Paleoceanography, 27, PA3228, doi:10.1029/2012PA002307.

M.O. Patterson, R. McKay, T. Naish, C. Escutia, F.J. Jimenez-Espejo, M.E. Raymo, M.E., 
S.R. Meyers, L. Tauxe, H. Brinkhuis, and IODP Expedition 318 Scientists,2014, 
\emph{Response of the East Antarctic Ice Sheet to orbital forcing during the Pliocene and Early Pleistocene}, Nature Geoscience, v. 7, p. 841-847.

D.J. Thomson, 1982, 
\emph{Spectrum estimation and harmonic analysis}: IEEE Proceedings, v. 70, p. 1055-1096.
}

\seealso{ \code{\link{confAdjust}},\code{\link{multiTest}},\code{\link{lowspec}}, \code{\link{mtm}}, \code{\link{mtmML96}}, \code{\link{mtmPL}}, and \code{\link{periodogram}} }


\examples{
\donttest{
# evaluate false positive rate for MTM-AR1 using AR1 surrogates
testBackground(noiseType="ar1",method="mtmAR1")

# evaluate false positive rate for MTM-AR1 using Power Law surrogates
testBackground(noiseType="pwrLaw",method="mtmAR1")
}
}