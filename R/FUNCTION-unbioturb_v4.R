### This function is a component of astrochron: An R Package for Astrochronology
### Copyright (C) 2022 Huaran Liu and Stephen Meyers

#'  Bioturbation removal function 
#'  
#'  unbioturb is a function to remove the boioturbation effects from a time series 
#'  given the bioturbation characteristics. It is an inverse process of function bioturb.
#'  Both processes model the bioturbation as a diffusive procedure (Guinasso and Schink 1975). 
#'  In unbioturb, the proxy series is deconvolved from an impulse response function determined
#'  by bioturbation characteristics (G = D/ML/v), building on the approach of Goreau (1980).
#'
#' @param G Control parameter in Guinasso and Schinck 1975. G = D/ML/v
#' @param ML Mix layer depth cm
#' @param v Sedimentation rate in cm/kyr
#' @param dat Input proxy time series. The order of dat follow the that of the depth: dat[1] corresponds to youngest point, dat[nt] corresponds to oldest point.
#' @param dtime Time resolution of dat in kyr
#' @param fhigh cut-off frequency for low pass filtering
#' @return return_list$dat_out  Bioturbed time series, same length as dat 
#' @return return_list$fc  Impulse response function

unbioturb <- function( dat, G, ML, v, pt = 0.2, wiener = TRUE, fhigh=NULL, output=1, genplot = TRUE, check = TRUE, verbose = TRUE)
{

   dat <- data.frame(dat)
   npts <- length(dat[,1])
   dtime <- dat[2,1]-dat[1,1]
   total_time <- (npts-1)*dtime

# error checking 
   if(check)
   {
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
  
  # ---- FFT deconvolution ----
  tdepo <- 500
  dat_out <- deconv(dat_interp$y-mean(dat_interp$y, na.rm = TRUE), fc/2, tdepo, dt, pt, wiener = wiener)+mean(dat_interp$y, na.rm = TRUE)
  
  # resample on the original sampling grid   
  out=approx(xout, dat_out, dat[,1], method = "linear")
  out=na.omit(as.data.frame(out))

# low pass filter for stabilization, if needed
  if(!is.null(fhigh)) 
   {
     if(fhigh<=(1/(2*dtime)) && fhigh >= (1/(length(out[,1])*dt))) out=taner(out, fhigh = fhigh, roll=10^10, check= FALSE, genplot = FALSE, verbose = FALSE)
     if(fhigh>(1/(2*dtime)) && verbose) cat("\n**** WARNING: fhigh set above Nyquist frequency.\n")
     if(fhigh<(1/(length(out[,1])*dtime)) && verbose) cat("\n**** WARNING: fhigh set below Rayleigh frequency.\n")
   }
   
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
    
    plot(out[,1],out[,2], xlim=c(min(dat[,1]),max(dat[,1])),ylim=c(min(dat[,2],out[,2]),max(dat[,2],out[,2])),
         main = "Unbioturbated Series (blue) & Proxy Observation (gray)", xlab = "Time (kyr, young -> old)",
         panel.first = grid(), type="l", col = "blue", lwd = 2)
    lines(dat[,1], dat[,2], type = "l", col = "#0000005F", lwd = 2)
    abline(v=dat_interp$x[length(dat_interp$x)*pt/2], col="black", lty=2, lwd=1)
    abline(v=dat_interp$x[length(dat_interp$x)*(1-pt/2)], col="black", lty=2, lwd=1)
  }
#  return_list <- list("fc" = fc, "dat_out" = dat_out)
#  return(return_list)
#  remove the portion of the series that is tapered, prior to output
  if (output==1) 
    {
      out <- subset(out, (out[1] > dat_interp$x[length(dat_interp$x)*pt/2]) & (out[1] < dat_interp$x[length(dat_interp$x)*(1-pt/2)]))
      if(verbose) 
        {
          cat(" * Tapered portion of series has been removed. \n")
          cat(" * New number of data points in unbioturbated series:",length(out[,1]),"\n")
        }  
      return(out)
    }     
  if (output==2) return(data.frame(cbind(xx,fc)))

}