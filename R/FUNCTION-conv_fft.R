### This function is a component of astrochron: An R Package for Astrochronology
### Copyright (C) 2021 Huaran Liu and Stephen Meyers

#' Convolution through Fast Fouriem Transform
#' 
#' x and y do not need to be of the same length. 
#'
#' @param x Input signal that needs to be convolved.
#' @param y Green function/ Impulse response function. 
#' @param dt time resolution for the input series
#' @param index index in the impulse response function that corresponds to the deposition time 
#' @return z Convolved output series. length(z) = length(x)
#' @import stats
#' @export

conv_fft <- function(x, y, index, dt){
  # Pad input series with 0 to the length of length(x)+length(y)-1
  nx <- length(x)
  ny <- length(y)
  nfft <- nx+ny-1
  
  x <- append(x, rep(0, ny-1))
  y <- append(y, rep(0, nx-1))
  
  fft_z <- Re(fft( fft(x)*fft(y), inverse = TRUE)/nfft)
  z <- fft_z[index: (nfft-(ny-index))]*dt
  return(z)
  
}