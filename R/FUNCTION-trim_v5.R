### This function is a component of astrochron: An R Package for Astrochronology
### Copyright (C) 2016 Stephen R. Meyers
###
###########################################################################
### trim: function to remove outliers - (SRM: March 9-10, 2012; June 30, 2012; 
###                                      May 20, 2013; August 6, 2015; July 22, 2016)
###
### automatically remove outliers from time series.
### outliers are determined using boxplot
###########################################################################
 
trim <- function (dat,c=1.5,genplot=T,verbose=T)
{

  if(verbose) cat("\n----- REMOVING OUTLIERS FROM STRATIGRAPHIC SERIES -----\n")
  
  ipts <- length(dat[,1]) 
  if(verbose) cat(" * Number of original data points=", ipts,"\n")

  eval <- !dat[,2] %in% boxplot.stats(dat[,2],coef=c)$out
  out <- rep(NA,ipts*2)
  dim(out) <- c(ipts,2)
  trimmed <- rep(NA,ipts*2)
  dim(trimmed) <- c(ipts,2)
  trimpts=0
  for (i in 1:ipts)
    { 
      if (isTRUE(eval[i])) 
        {
          out[i,1]= dat[i,1]
          out[i,2]=dat[i,2]
          trimpts = trimpts + 1
        }   
      else 
        {
          trimmed[i,1]= dat[i,1]
          trimmed[i,2]= dat[i,2]
        }   
     }
  
    if(verbose) cat(" * Number of data points post-trimming=", trimpts,"\n")
### grab all data except 'NA' values, assign to dat2
    dat2 <- data.frame(subset(out,!(out[,2] == 'NA')))

### grab all trimmed data except 'NA' values, assign to datrimmed
    datrimmed <- data.frame(subset(trimmed,!(trimmed[,2] == 'NA')))

    if(genplot)
     {
### plot data series. Note, cex is the factor by which to increase or decrease default symbol size
      par(mfrow=c(2,2))
      xrange = range(c(dat[,1]))
      yrange = range(c(dat[,2]))
      plot(dat2, cex=.5,xlim=xrange,ylim=yrange,xlab="Location",ylab="Value",main="Data Series")
      lines(dat2)
      par(new=TRUE)
      plot(datrimmed, cex=.5, xlim=xrange,ylim=yrange,col="red")  
     
### plot the denisty and the histogram together
      hist(dat2[,2],freq=F,xlab="Value",main="Distribution values"); lines(density(dat2[,2], bw="nrd0"),col="red"); grid()
### boxplot
      boxplot(dat2[,2])
### Normal probabilty plot (Normal Q-Q Plot)
      qqnorm(dat2[,2]); qqline(dat2[,2], col="red"); grid()
    }

    return(dat2)
    
### END function trim
}
