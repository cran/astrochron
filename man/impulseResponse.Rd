\name{impulseResponse}
\alias{impulseResponse}
\title{Impulse response function calculation}
\description{
Calculate the analytical response function from an impulse forcing using the 1-D advection diffusion model proposed 
in Schink and Guinasso (1975) to model the bioturbation impact on climate time series.
}
\usage{
impulseResponse(G, ML = NULL, v = NULL, nt = 500, genplot = FALSE,
  verbose = FALSE)
}
\arguments{
  \item{G}{Bioturbation parameter. G = D/ML/v}
  \item{ML}{Mix layer depth (cm)}
  \item{v}{Sedimentation rate (cm/kyr)}
  \item{nt}{Number of steps after the signal is deposited.}
  \item{genplot}{Generate summary plots? (T or F)}
  \item{verbose}{Verbose output? (T or F)}
}
\value{
fc Impulse response function
}
\references{
Guinasso, N.L. and Schinck, D.R., 1975, \emph{Quantitative estimates of biological mixing rates in abyssal sediments}, J. Geophys. Res., 80, 3032-3043.

Goreau, T.J., 1977, \emph{Quantitative effects of sediment mixing on stratigraphy and biogeochemistry: a signal theory approach}, Nature, 256, 730-732.

Goreau, T.J., 1980, \emph{Frequency sensitivity of the deep-sea climatic record}, Nature, 287, 620-622.

Liu, H., Meyers, S.R., and Marcott, S.A., 2021, \emph{Unmixing dee-sea paleoclimate records: A study on bioturbation effects through convolution and deconvolution}, Earth and Planetary Science Letters, 564, 116883.
}
\examples{
G <- 4
ML <- 10
v <- 1
# take a look at the IRF 
impulseResponse(G=4, ML = 10, v = 1, genplot = TRUE)
}