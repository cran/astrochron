### This function is a component of astrochron: An R Package for Astrochronology
### Copyright (C) 2021 Huaran Liu and Stephen Meyers

#' Wiener Deconvolution through Fast Fouriem Transform. cosin taper is applied to remove edge effects
#' signal-to-noise ratio  is chosen to be 0.05, gamma is also chosen to be 0.05
#' x and y do not need to be of the same length. 
#'
#' @param x Bioturbed time series that needs to be deconvolved.
#' @param y Green function/ Impulse response function. 
#' @param dt time resolution for the input series
#' @param index index in the impulse response function that corresponds to the deposition time 
#' @return z deconvolved/un-bioturbed  series. length(z) = length(x)
#' @import stats astrochron
#' @export

deconv <- function(x, y, index, dt, pt = 0.2, wiener = TRUE){
  # Pad input series with 0 to the length of length(x)+length(y)-1
  nx <- length(x)
  ny <- length(y)
  nfft <- nx+ny-1
  
  # apply cosTaper to x before being padded, note cosTaper takes data frame as input
  x.df <- data.frame(1:length(x), x)
  x <- cosTaper(x.df, p = pt, genplot = F, verbose = F)$x
  x <- c( rep(0, index-1),x,rep(0, ny-index) )/dt
  y <- append(y, rep(0, nx-1))
  
 
  y.fft <- fft(y)
  y.fft.mod <- abs(y.fft)^2

  if(wiener){
  # Wiener decovolution algorithm
  snr <- rep(1, length(y.fft))
  snr[y.fft.mod<0.05] <- 0.05
  fft_z <- Re(fft(    fft(x)/fft(y)*abs(fft(y))^2/(abs(fft(y))^2+1/snr), inverse = TRUE)/nfft)
  }else{
  fft_z <- Re(fft(    fft(x)/fft(y), inverse = TRUE)/nfft)
  }
  
  z <- fft_z[1:nx]
  return(z)
  
}