\name{multiTest}
\alias{multiTest}
\title{Adjust spectral p-values for multiple comparisons}
\description{
Adjust spectral p-values for multiple comparisons, using a range of approaches.
}
\usage{
multiTest(spec,flow=NULL,fhigh=NULL,pl=T,output=T,genplot=T,verbose=T)
}
\arguments{
  \item{spec}{A data frame with two columns: frequency, uncorrected confidence level. If more than 2 columns are input, the results are assumed to come from periodogram, mtm, mtmML96, lowspec or mtmPL.}
  \item{flow}{Vector of lower bounds for each frequency band of interest. Order must match fhigh.}
  \item{fhigh}{Vector of upper bounds for each frequency band of interest. Order must match flow.}
  \item{pl}{Include graphs of uncorrected p-values? (T or F)}
  \item{output}{Return results as new data frame? (T or F)}
  \item{genplot}{Generate summary plots? (T or F)}
  \item{verbose}{Verbose output? (T or F)}
}
\details{
Multiple testing is a common problem in the evaluation of power spectrum peaks 
(Vaughan et al., 2011; Crampton et al., PNAS).  To address the issue
of multiple testing, a range of approaches have been advocated. This function will conduct an 
assessment using six approaches: Bonferroni, Holm (1979), Hochberg (1998), Hommel (1988),
Benjamini & Hochberg (1995, a.k.a. "false discovery rate"), 
Benjamini & Yekutieli (2001, a.k.a. "false discovery rate").  See the function p.adjust
for additional information on these six approaches.

In conducting these assessments, it is important that the spectral analysis is conducted 
without zero-padding. If one is exclusively concerned with particular frequency bands a priori 
(e.g., those associated with Milankovitch cycles), the statistical power of the method
can be improved by restricting the analysis to those frequency bands (use options 'flow' and 'fhigh').

Application of these multiple testing corrections does not guarantee that the spectral background
is appropriate.  To address this issue, carefully examine the fit of the spectral background,
and also conduct simulations with the function testBackground.
}

\references{
Y. Benjamini, and Y. Hochberg, 1995,
 \emph{Controlling the false discovery rate: a practical and powerful approach to multiple testing}. Journal of the Royal Statistical Society Series B 57, 289-300.

Y. Benjamini, and D. Yekutieli, 2001, 
 \emph{The control of the false discovery rate in multiple testing under dependency}. Annals of Statistics 29, 1165-1188.

J.S. Campton, S.R. Meyers, R.A. Cooper, P.M Sadler, M. Foote, D. Harte, 2018, 
\emph{Pacing of Paleozoic macroevolutionary rates by Milankovitch grand cycles}: Proceedings of the National Academy of Sciences, doi:10.1073/pnas.1714342115.

S. Holm, 1979, 
 \emph{A simple sequentially rejective multiple test procedure}. Scandinavian Journal of Statistics 6, 65-70.

G. Hommel, 1988, 
 \emph{A stagewise rejective multiple test procedure based on a modified Bonferroni test}. Biometrika 75, 383-386.

Y. Hochberg, 1988, 
 \emph{A sharper Bonferroni procedure for multiple tests of significance}. Biometrika 75, 800-803.

J.P. Shaffer, 1995,
 \emph{Multiple hypothesis testing. Annual Review of Psychology 46, 561-576}. (An excellent review of the area.)

S. Vaughan, R.J. Bailey, and D.G. Smith, 2011, 
 \emph{Detecting cycles in stratigraphic data: Spectral analysis in the presence of red noise}. Paleoceanography 26, PA4211, doi:10.1029/2011PA002195.
}

\seealso{ \code{\link{p.adjust}},\code{\link{testBackground}},\code{\link{confAdjust}},\code{\link{lowspec}}, \code{\link{mtm}}, \code{\link{mtmML96}}, \code{\link{mtmPL}}, and \code{\link{periodogram}} }


\examples{
# generate example series with periods of 400 ka, 100 ka, 40 ka and 20 ka
ex = cycles(freqs=c(1/400,1/100,1/40,1/20),start=1,end=1000,dt=5)

# add AR1 noise
noise = ar1(npts=200,dt=5,sd=.5)
ex[2] = ex[2] + noise[2]

# first, let's look at mtm with conventional AR1 background
spec=mtm(ex,padfac=1,ar1=TRUE,output=1)

# when blindly prospecting for cycles, it is necessary to consider all of the 
#  observed frequencies in the test
multiTest(cb(spec,c(1,4)),output=FALSE)

# if, a priori, you are only concerned with the Milankovitch frequency bands, 
#  restrict your analysis to those bands (as constrained by available sedimentation
#  rate estimates and the frequency resolution of the spectrum). in the example below, 
#  the mtm bandwidth resolution is employed to search frequencies nearby the 
#  Milankovitch-target periods.
flow=c((1/400)-0.003,(1/100)-0.003,(1/41)-0.003,(1/20)-0.003)
fhigh=c((1/400)+0.003,(1/100)+0.003,(1/41)+0.003,(1/20)+0.003)
multiTest(cb(spec,c(1,4)),flow=flow,fhigh=fhigh,output=FALSE)

# now try with the lowspec method. this uses prewhitening, so it has one less data point.
spec=lowspec(ex,padfac=1,output=1)
flow=c((1/400)-0.003015075,(1/100)-0.003015075,(1/41)-0.003015075,(1/20)-0.003015075)
fhigh=c((1/400)+0.003015075,(1/100)+0.003015075,(1/41)+0.003015075,(1/20)+0.003015075)
multiTest(cb(spec,c(1,4)),flow=flow,fhigh=fhigh,output=FALSE)

# for comparison...
multiTest(cb(spec,c(1,4)),output=FALSE)
}