### This function is a component of astrochron: An R Package for Astrochronology
### Copyright (C) 2014 Stephen R. Meyers
###
###########################################################################
### xplot: create cross plot of two variable, with density estimates along
###        axes - (SRM: April 16, 2014)
###
###########################################################################


xplot <- function (x,y,xlab=NULL,ylab=NULL,main=NULL)
{

  cat("\n----- GENERATING CROSS-PLOT WITH DENSITY ESTIMATES ON AXES -----\n")
  x<-data.frame(x)
  y<-data.frame(y)
  if( length(x[,1]) != length(y[,1]) ) stop("***** ERROR: number of data points in two series is not equivalent!")

  xmax=max(x[,1])
  xmin=min(x[,1])
  ymax=max(y[,1])
  ymin=min(y[,1])

  par(fig=c(0,0.8,0,0.8))
  plot(x[,1],y[,1],xlim=c(xmin,xmax),ylim=c(ymin,ymax),cex=1.2,cex.axis=1.2,cex.lab=1.2,xlab="",ylab="")
  if(is.null(xlab)) xlab=colnames(x)
  if(is.null(ylab)) ylab=colnames(y)
  if(is.null(main)) main=""
  mtext(ylab,side=2,line=2.5,cex=1.2)
  mtext(xlab,side=1,line=2.5,cex=1.2)
  par(fig=c(0,0.8,0.55,1), new=TRUE)
  z <- density(x[,1], bw="nrd0")
  plot(z$x,z$y, type='l',axes=FALSE,xlab="",ylab="",lwd=2,xlim=c(xmin,xmax))
  par(fig=c(0.65,1,0,0.8),new=TRUE)
  z <- density(y[,1], bw="nrd0")
  plot(z$y, z$x, type='l',axes=FALSE,xlab="",ylab="",lwd=2,ylim=c(ymin,ymax))
  mtext(main, side=3, cex=1.5, outer=TRUE, line=-3)

### END function xplot
}
