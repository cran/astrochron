\name{deconv}
\alias{deconv}
\title{Wiener Deconvolution through Fast Fourier Transform}
\description{
Wiener Deconvolution through Fast Fourier Transform. 
}
\usage{
deconv(x, y, index, dt, pt = 0.2, wiener = TRUE)
}
\arguments{
  \item{x}{Time series that needs to be deconvolved.}
  \item{y}{Green function/Impulse response function.}
  \item{index}{index in the impulse response function that corresponds to the deposition time.}
  \item{dt}{time resolution for the input series.}
  \item{pt}{Cosine-tapered window parameter: pt is the percent of the data series tapered (choose 0-1). When pt=1, this is equivalent to a Hann taper.}
  \item{wiener}{Apply Wiener filter? (T or F)}
}
\value{
z deconvolved/un-bioturbed  series. length(z) = length(x)
}
\details{
function 'deconv' is used by function 'unbioturb'. A cosine taper is applied to remove edge effects.
The signal-to-noise ratio is chosen to be 0.05, and gamma is also chosen to be 0.05.
x and y do not need to be of the same length. For additional information see Liu et al. (2021)
}
\references{
Liu, H., Meyers, S.R., and Marcott, S.A., 2021, \emph{Unmixing dee-sea paleoclimate records: A study on bioturbation effects through convolution and deconvolution}, Earth and Planetary Science Letters, 564, 116883.
}
