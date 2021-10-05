\name{bioturb}
\alias{bioturb}
\title{Bioturbate time series using diffusion model from Guinasso and Schinck (1975), as in Liu et al. (2021)}
\description{
'bioturb' is a function for simulating the bioturbed time series when bioturbation is modeled as a diffusive process.
It implements the methodology outlined in Liu et al. (2021), which builds on the approaches of (Guinasso and Schink, 1975),
Goreau (1977), and Goreau (1980). Given the bioturbation parameters, an impulse response function is first calculated 
using function 'impulseResponse'. This function is then convolved with the input time series (true signal) to output the 
bioturbed time series. Note that the input true signal 'dat' and impulse response function are interpolated to the same 
resolution before convolution.
}
\usage{
bioturb(dat, G, ML, v, output=1, genplot=TRUE, verbose=TRUE)
}
\arguments{
  \item{dat}{Stratigraphic series to be bioturbated. First column should be age (kyr), second column should be data value.}
  \item{G}{Control parameter in Guinasso and Schinck, 1975. G = D/ML/v}
  \item{ML}{Mix layer depth (cm)}
  \item{v}{Sedimentation rate (cm/kyr)}
  \item{output}{Which results would you like to return to the console? (0) no output; (1) return bioturbated series; (2) return impulse response}
  \item{genplot}{Generate summary plots? (T or F)}
  \item{verbose}{Verbose output? (T or F)}
}

\references{
Guinasso, N.L. and Schinck, D.R., 1975, \emph{Quantitative estimates of biological mixing rates in abyssal sediments}, J. Geophys. Res., 80, 3032-3043.

Goreau, T.J., 1977, \emph{Quantitative effects of sediment mixing on stratigraphy and biogeochemistry: a signal theory approach}, Nature, 256, 730-732.

Goreau, T.J., 1980, \emph{Frequency sensitivity of the deep-sea climatic record}, Nature, 287, 620-622.

Liu, H., Meyers, S.R., and Marcott, S.A., 2021, \emph{Unmixing dee-sea paleoclimate records: A study on bioturbation effects through convolution and deconvolution}, Earth and Planetary Science Letters, 564, 116883.
}

\examples{
# as a test series, use the three dominant precession terms from Berger et al. (1992)
ex1=cycles()

# mix it
res1 <- bioturb(ex1, G=4, ML=10, v=1, genplot = TRUE)

# un-mix it
res2=unbioturb(res1, G=4, ML=10, v=1, genplot = TRUE)

pl()
plot(ex1,type="l",main="black=signal, blue=bioturbated, red=unbioturbated",lwd=3)
lines(res2,col="red")
lines(res1,col="blue")
}