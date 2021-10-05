\name{conv_fft}
\alias{conv_fft}
\title{Convolution through Fast Fourier Transform}
\description{
Convolution through Fast Fourier Transform
}
\usage{
conv_fft(x, y, index, dt)
}
\arguments{
  \item{x}{Input signal that needs to be convolved}
  \item{y}{Green's function; impulse response function}
  \item{index}{index in the impulse response function that corresponds to the deposition time}
  \item{dt}{time resolution for the input series}
}
\value{
z Convolved output series. length(z) = length(x)
}
\details{
Function 'conv_fft' is used by function 'bioturb'. x and y do not need to be of the same length.
}