\name{confAdjust}
\alias{confAdjust}
\title{Adjust spectrum confidence levels for multiple comparisons}
\description{
Adjust spectrum confidence levels for multiple comparisons, using the Bonferroni correction
}
\usage{
confAdjust(spec,npts,dt,tbw=3,ntap=5,flow=NULL,fhigh=NULL,output=T,
    xmin=df,xmax=NULL,pl=1,genplot=T,verbose=T)
}
\arguments{
  \item{spec}{A data frame with three columns: frequency, power, background power. If 8 columns are input, the results are assumed to come from mtm, mtmML96, lowspec or mtmPL. If 9 columns are input, the results are assumed to come from periodogram.}
  \item{npts}{Number of points in stratigraphic series.}
  \item{dt}{Sampling interval of stratigraphic series.}
  \item{tbw}{MTM time-bandwidth product.}
  \item{ntap}{Number of DPSS tapers to use.}
  \item{flow}{Vector of lower bounds for each frequency band of interest. Order must match fhigh.}
  \item{fhigh}{Vector of upper bounds for each frequency band of interest. Order must match flow.}
  \item{output}{Output data frame? (T or F)}
  \item{xmin}{Smallest frequency for plotting.}
  \item{xmax}{Largest frequency for plotting.}
  \item{pl}{Plotting option (1-4): 1=linear frequency & log power, 2=log frequency & power, 3=linear frequency & power, 4=log frequency & linear power.}
  \item{genplot}{Generate summary plots? (T or F)}
  \item{verbose}{Verbose output? (T or F)}
}
\details{
Multiple testing is a common problem in the evaluation of power spectrum peaks 
(Vaughan et al., 2011; Crampton et al., PNAS).  To address the issue of multiple testing, 
a range of approaches have been advocated. This function will conduct an assessment 
using the Bonferroni correction, which is the simplest, and also the most conservative,
of the common approaches (it is overly pessimistic).

If one is exclusively concerned with particular frequency bands a priori 
(e.g., those associated with Milankovitch cycles), the statistical power of the method
can be improved by restricting the analysis to those frequency bands (use options 'flow' and 'fhigh').

Application of multiple testing corrections does not guarantee that the spectral background
is appropriate.  To address this issue, carefully examine the fit of the spectral background,
and also conduct simulations with the function testBackground.
}
\references{
J.S. Campton, S.R. Meyers, R.A. Cooper, P.M Sadler, M. Foote, D. Harte, 2018, 
\emph{Pacing of Paleozoic macroevolutionary rates by Milankovitch grand cycles}: Proceedings of the National Academy of Sciences, doi:10.1073/pnas.1714342115.

S. Vaughan, R.J. Bailey, and D.G. Smith, 2011, 
 \emph{Detecting cycles in stratigraphic data: Spectral analysis in the presence of red noise}. Paleoceanography 26, PA4211, doi:10.1029/2011PA002195.
}

\seealso{ \code{\link{testBackground}},\code{\link{multiTest}},\code{\link{spec.mtm}}, \code{\link{lowspec}}, and \code{\link{periodogram}} }

\examples{
# generate example series with periods of 400 ka, 100 ka, 40 ka and 20 ka
ex = cycles(freqs=c(1/400,1/100,1/40,1/20),start=1,end=1000,dt=5)

# add AR1 noise
noise = ar1(npts=200,dt=5,sd=.5)
ex[2] = ex[2] + noise[2]

# first, let's use mtm with conventional AR1 background
spec=mtm(ex,padfac=1,ar1=TRUE,output=1)

# when blindly prospecting for cycles, it is necessary to consider all of the 
#  observed frequencies in the test
confAdjust(spec,npts=200,dt=5,tbw=3,ntap=5,output=FALSE)

# if, a priori, you are only concerned with the Milankovitch frequency bands, 
#  restrict your analysis to those bands (as constrained by available sedimentation
#  rate estimates and the frequency resolution of the spectrum). in the example below, 
#  the mtm bandwidth resolution is employed to search frequencies nearby the 
#  Milankovitch-target periods.
flow=c((1/400)-0.003,(1/100)-0.003,(1/41)-0.003,(1/20)-0.003)
fhigh=c((1/400)+0.003,(1/100)+0.003,(1/41)+0.003,(1/20)+0.003)
confAdjust(spec,npts=200,dt=5,tbw=3,ntap=5,flow=flow,fhigh=fhigh,output=FALSE)

# now try with the lowspec method. this uses prewhitening, so it has one less data point.
spec=lowspec(ex,padfac=1,output=1)
flow=c((1/400)-0.003015075,(1/100)-0.003015075,(1/41)-0.003015075,(1/20)-0.003015075)
fhigh=c((1/400)+0.003015075,(1/100)+0.003015075,(1/41)+0.003015075,(1/20)+0.003015075)
confAdjust(spec,npts=199,dt=5,tbw=3,ntap=5,flow=flow,fhigh=fhigh,output=FALSE)

# for comparison...
confAdjust(spec,npts=199,dt=5,tbw=3,ntap=5,output=FALSE)
}