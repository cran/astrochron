\name{unbioturb}
\alias{unbioturb}
\title{
Bioturbation removal function following the approach of Liu et al (2021)
} 
\description{
 'unbioturb' is a function to remove bioturbation effects from a time series 
 given the bioturbation parameters. It implements the method outlined in 
 Liu et al. (2021), which builds on the approaches of (Guinasso and Schink, 1975),
 Goreau (1977), and Goreau (1980). 'unbioturb' is the inverse of the function 'bioturb', 
 both of which model bioturbation as a diffusive process (Guinasso and Schink, 1975). 
 In 'unbioturb', the proxy series is deconvolved from an impulse response function 
 determined by the bioturbation characteristics, G = D/ML/v.
}
\usage{
unbioturb(dat, G, ML, v, pt = 0.2, wiener = TRUE, fhigh=NULL, output = 1,
  genplot = TRUE, check = TRUE, verbose = TRUE)
}
\arguments{
  \item{dat}{Stratigraphic series to be bioturbated. First column should be age (kyr), second column should be data value.}
  \item{G}{Control parameter in Guinasso and Schinck, 1975. G = D/ML/v}
  \item{ML}{Mix layer depth (cm)}
  \item{v}{Sedimentation rate (cm/kyr)}
  \item{pt}{Cosine-tapered window parameter: pt is the percent of the data series tapered (choose 0-1). When pt=1, this is equivalent to a Hann taper.}
  \item{wiener}{Apply Wiener filter for deconvolution stabilization? (T or F)}
  \item{fhigh}{Taner filter cut-off frequency for deconvolution stabilization. By default, no Taner lowpass filter is applied.}
  \item{output}{Which results would you like to return to the console? (0) no output; (1) return unbioturbated series; (2) return impulse response}
  \item{genplot}{Generate summary plots? (T or F)}
  \item{check}{Conduct compliance checks before processing? (T or F) In general this should be activated; the option is included for Monte Carlo simulation.}
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
