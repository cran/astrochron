### This function is a component of astrochron: An R Package for Astrochronology
### Copyright (C) 2021 Huaran Liu and Stephen Meyers

#' Impulse response function calculation
#' 
#' Calculate the analytical response function from an impulse forcing using the 1-D advection diffusion model proposed 
#' in Schink and Guinasso 1975 to model the bioturbation impact on climate time series
#' @param nt  Number of steps after the signal is deposited.
#' @param G Bioturbation parameter. G = D/ML/v
#' @return fc Impulse response function
#' @useDynLib TestBioturb Root_Search
#' @useDynLib TestBioturb Impulse_Response
#' @export

impulseResponse <- function(G, ML = NULL, v = NULL,  nt = 500, genplot = FALSE, verbose = FALSE){
  
  # Wrapper for root_search C function
  root_search <- function(G){
    .Call("Root_Search", G)
  }
  # Wrapper for Impulse_Response C function
  response <- function(dt, nt, G, a, verbose = FALSE){
    v <- 0
    if(verbose){
      v <- 1
    }
    .Call("Impulse_Response", dt, as.integer(nt), G, alpha, as.integer(v))
  }
  
  
  dt <- 0.01
  Nroots <- 300
  alpha <- root_search(G)
  fc <- response(dt, nt, G, alpha, verbose)
  
  # ploting
  if(genplot){
    if (is.null(ML) | is.null(v)){
      
      xx = seq(dt, 10, dt)
      titleA = sprintf("Impulse Responsen Series ( G = %3.2f)", G)
      plot(xx, fc, type = "l", cex = 0.5, lwd = 2,
         xlab="Depth/Time",ylab="Value", main=titleA,
         panel.first = grid())
      abline(v=500*dt, col="red", lty=2, lwd=3)}
    else{
      par(mfrow=c(3,1))
      
      xx = seq(dt, 10, dt)
      titleA = sprintf("G = %3.2f, ML = %3.2f cm, v = %3.2f cm/kyr", G, ML, v)
      plot(xx, fc, type = "l", cex = 0.5, lwd = 2,
           xlab="Depth/Time (1 unit coresponds to mix layer depth)",
           ylab="Value",main = paste("Impulse Response Series", "\n", titleA),
           panel.first = grid())
      abline(v=500*dt, col="red", lty=2, lwd=3)
      
      
      xx = seq(dt, 10, dt)*ML
      plot(xx, fc, type = "l", cex = 0.5, lwd = 2,
           xlab="Depth (cm)",ylab="Value", 
           main = paste("Impulse Response Series", "\n", titleA),
           panel.first = grid())
      abline(v=500*dt*ML, col="red", lty=2, lwd=3)
      
      xx = seq(dt, 10, dt)*ML/v
      plot(xx, fc, type = "l", cex = 0.5, lwd = 2,
           xlab="Time (kyr,  Present    --->       Past) ",ylab="Value",
           main= paste("Impulse Response Series", "\n", titleA),
           panel.first = grid())
      abline(v=500*dt*ML/v, col="red", lty=2, lwd=3)
      
    }
  }
  return(fc)
}