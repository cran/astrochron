### This function is a component of astrochron: An R Package for Astrochronology
### Copyright (C) 2019 Stephen R. Meyers
###
###########################################################################
### repl0: replace negative data values with 0 - (SRM: May 2, 2012; April 30, 2013; 
###                                               May 20, 2013; July 22, 2016; 
###                                               June 19, 2019)
###
###########################################################################

repl0 <- function (dat, ID=T, genplot=T, verbose=T)
{

  if(verbose) cat("\n----- REPLACING NEGATIVE DATA VALUES WITH 0 -----\n")

  dat <- data.frame(dat)
  ipts <- length(dat[,1]) 
  if(verbose) cat(" * Number of data points=", ipts,"\n")
  cols=length(dat)
  if(verbose) cat(" * Number of columns=", cols,"\n")
 
  if( any(is.na(dat[,2:cols])) ) cat("**** WARNING: empty entries, NA or NaN are present in data frame!\n")
  numneg = sum(dat[,2:cols]<0, na.rm=TRUE)
  if(verbose) cat(" * Number of negative values=", numneg,"\n")
  
  repl.0 <- function(x) {x[x<0] <- 0; x}
  dat[,2:cols] <- repl.0(dat[,2:cols])
  
 if(cols == 2 && genplot)
  {
### plots
   par(mfrow=c(2,2))
   plot(dat,cex=0.5,xlab="Location",ylab="Value",main="Data Series"); lines(dat)
### plot the denisty and the histogram together
   hist(dat[,2],freq=F,xlab="Value",main="Distribution of isoloated values") 
# skip density plot if NaN present in data frame
   if( !any(is.na(dat[,2])) ) lines(density(dat[,2], bw="nrd0"),col="red")
### boxplot
   boxplot(dat[,2],ylab="Value",main="Boxplot for isolated values")
### Normal probabilty plot (Normal Q-Q Plot)
   qqnorm(dat[,2]); qqline(dat[,2], col="red")
  }
  
  return(dat)

### END function repl0
}
