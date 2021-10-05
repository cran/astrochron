### This function is a component of astrochron: An R Package for Astrochronology
### Copyright (C) 2021 Huaran Liu and Stephen Meyers

#' Bioturbed time series using diffusion model from Guinasso and Schinck 1975
#' 
#' bioturb is a function for simulating a bioturbed time series when bioturbation is modeled as a diffusion process.
#' Given the bioturbation parameters, an impulse response function is first calculated using function impulseResponse. 
#' This function is then convolved with the input time series (true signal) to output the bioturbed time series.
#' Note that the input true signal dat and impulse response function are interpolated to the same resolution before 
#' convolution.   
#' @param G Control parameter in Guinasso and Schinck 1975. G = D/ML/v
#' @param ML Mix layer depth cm
#' @param v Sedimentation rate in cm/kyr
#' @param dat Input time series
#' @param total_time Total length of dat in kyr
#' @return return_list$dat_out  Bioturbed time series, same length as dat 
#' @return return_list$fc  Impulse response function#' 


bioturb <- function( dat, G, ML, v, output=1, genplot = TRUE, verbose=TRUE)
{

   dat <- data.frame(dat)
   npts <- length(dat[,1])
   dtime <- dat[2,1]-dat[1,1]
   total_time <- (npts-1)*dtime

# error checking 
   if(dtime<0)
     { 
       if (verbose) cat("\n * Sorting data into increasing height/depth/time, removing empty entries\n")
       dat <- dat[order(dat[,1], na.last = NA, decreasing = F), ]
       dtime <- dat[2,1]-dat[1,1]
       npts <- length(dat[,1])
     }

   dtest <- dat[2:npts,1]-dat[1:(npts-1),1] 
   epsm=1e-9
   if( (max(dtest)-min(dtest)) > epsm ) 
     {
       cat("\n**** ERROR: sampling interval is not uniform.\n")
       stop("**** TERMINATING NOW!")
     }

   if (verbose) 
    {
      cat(" * Number of data points in stratigraphic series:",npts,"\n")
      cat(" * Stratigraphic series length (space or time):",total_time,"\n")
      cat(" * Sampling interval (space or time):",dtime,"\n")
    }

    
  # Get impulse response function, interpolate and then convolve
  fc <- impulseResponse(G)
  # impulse reponse is at 0.01 (normalized) resolution, interpolate data so it is on same resolution grid
  dt <- 0.01
  dtime0 <- ML/v*dt
  xout <- seq(dat[1,1], dat[npts,1], dtime0)
  
  dat_interp <- approx(dat[,1], dat[,2], xout, method = "linear", rule = 2)
  
  tdepo <- 500  
  dat_out <- conv_fft(dat_interp$y, fc/2, tdepo, dt)

  # resample on the original sampling grid   
  out=approx(xout, dat_out, dat[,1], method = "linear")
  out=na.omit(as.data.frame(out))
  
  if(genplot){
    layout(matrix(c(1,2,3,3), 2, 2, byrow = TRUE))
    
    titleA = sprintf("G = %3.2f, ML = %3.2f cm, v = %3.2f cm/kyr", G, ML, v)
    
    xx = seq(dt, 10, dt)
    plot(xx, fc, type = "l", cex = 0.5, lwd = 2,
         xlab="Depth/Time (1 unit coresponds to mix layer depth)",
         ylab="Value",main = paste("Impulse Response Series", "\n", titleA),
         panel.first = grid())
    abline(v=tdepo*dt, col="red", lty=2, lwd=3)
    
    xx = seq(dt, 10, dt)*ML/v
    plot(xx, fc, type = "l", cex = 0.5, lwd = 2,
         xlab="Time (kyr,  Present    --->       Past) ",ylab="Value",
         main= paste("Impulse Response Series", "\n", titleA),
         panel.first = grid())
    abline(v=tdepo*dt*ML/v, col="red", lty=2, lwd=3)
    
    plot(dat[,1], dat[,2], type = "l", col = "black", lwd = 2,
              main = "Input (black) & Output (blue)", xlab = "Time (kyr, young -> old)",
              panel.first = grid())
    lines(out[,1], out[,2], type = "l", col = "blue", lwd = 2)
  }
#  return_list <- list("fc" = fc, "dat_out" = dat_out)
#  return(return_list)
  if (output==1) return(out)
  if (output==2) return(data.frame(cbind(xx,fc)))

}