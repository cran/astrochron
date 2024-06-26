### This function is a component of astrochron: An R Package for Astrochronology
### Copyright (C) 2023 Stephen R. Meyers
###
###########################################################################
### idPts: identify points in plot (SRM: Oct. 29, 2012; Nov. 23, 2012;
###                                      May 20, 2013; June 10-11, 2013;
###                                      June 13, 2013; July 27, 2013;
###                                      April 10, 2014; April 21, 2014;
###                                      April 23, 2014; Nov. 19, 2014;
###                                      February 4, 2015; June 1, 2015;
###                                      September 14, 2017; September 12, 2023)
###########################################################################


idPts <- function (dat1,dat2=NULL,ptsize=0.5,xmin=NULL,xmax=NULL,ymin=NULL,ymax=NULL,logx=F,logy=F,plotype=1,annotate=1,iso=F,output=1,verbose=T)
{    
    if(verbose) cat("\n----- INTERACTIVELY IDENTIFY POINTS IN PLOT -----\n")
    dat1 <- data.frame(dat1)
    if(!is.null(dat2)) dat2 <- data.frame(dat2)
    if(length(dat1)==1 || length(dat1)==2) setIn=1
    if(length(dat1)>=3) setIn=2
    
    if(length(dat1)>3) cat("\n  WARNING:", length(dat1),"columns detected in dat1. Will only use first three columns.\n\n")

# added option to have three columns in dat1: location, variable x, variable y
    if(length(dat1) < 2 && is.null(dat2)) stop("**** TERMINATING: you must input 2 variables")
    if(length(dat1) < 2 && !is.null(dat2)) 
      {
        if(length(dat1) != 1 || length(dat2) != 1) stop("**** TERMINATING: you must input 2 variables")  
      }
    if(length(dat1) >= 2 && !is.null(dat2)) stop("**** TERMINATING: dat1 has >= 2 variables; do not use dat2")

    if(setIn == 1)
     {
      if(is.null(dat2))
        {
          xx <- dat1[,1]
          yy <- dat1[,2]
        } 
      if(!is.null(dat2))
        { 
          xx <- dat1[,1]
          yy <- dat2[,1]
        }  
     }

    if(setIn == 2)
      {
        loc = dat1[,1]
        xx <- dat1[,2]
        yy <- dat1[,3] 
      }     
   
    if(length(xx) < 2) stop("**** TERMINATING: input must have more than one point")
    if(length(xx) != length(yy)) stop("**** TERMINATING: the two variables must have the same number of points")

    if(verbose) cat(" * Select points by clicking\n")
    if(verbose) cat("   Stop by pressing ESC-key (Mac) or STOP button (Windows)\n")

    if(is.null(xmin)) xmin=min(xx)
    if(is.null(xmax)) xmax=max(xx)
    if(is.null(ymin)) ymin=min(subset(yy, (xx >= xmin) & (xx <= xmax) ) )
    if(is.null(ymax)) ymax=max(subset(yy, (xx >= xmin) & (xx <= xmax) ) )
   
    par(mfrow=c(1,1))
    if(logx && logy) setlog="xy"
    if(logx && !logy) setlog="x"
    if(!logx && logy) setlog="y"
    if(!logx && !logy) setlog=""
    if (plotype == 1) { plot(xx,yy, main="Select data points (press ESC-key or STOP to exit)",xlim=c(xmin,xmax),ylim=c(ymin,ymax),bty="n",lwd=2,cex.axis=1.1,cex.lab=1.1,cex=ptsize,xlab="",ylab="",log=setlog); lines(xx,yy,col="red") }
    if (plotype == 2) { plot(xx,yy, main="Select data points (press ESC-key or STOP to exit)",xlim=c(xmin,xmax),ylim=c(ymin,ymax),bty="n",cex.axis=1.1,cex.lab=1.1,cex=ptsize,log=setlog) }
    if (plotype == 3) { plot(xx,yy, type="l", main="Select data points (press ESC-key or STOP to exit)",xlim=c(xmin,xmax),ylim=c(ymin,ymax),bty="n",lwd=2,cex.axis=1.1,cex.lab=1.1,log=setlog) }

## this script modified from '?identify' in R
identifyPch <- function(x, y=NULL, loc=NULL, n=length(x), pch=19, cex, ...)
{
    xy <- xy.coords(x, y); x <- xy$x; y <- xy$y
    sel <- rep(FALSE, length(x)); res <- integer(0)
    while(sum(sel) < n) { 
        ans <- identify(x[!sel], y[!sel], n=1, plot=F, ...)
        if(!length(ans)) break
        ans <- which(!sel)[ans]
        points(x[ans], y[ans], pch = pch, cex = cex, col="blue")
        sel[ans] <- TRUE
        res <- c(res, ans)
        if(annotate==1) text(x[ans],y[ans],round(x[ans],digits=2),col = "blue", pos = 3, offset = 1.2, font=2, xpd = TRUE)
        if(annotate==1) text(x[ans],y[ans],round(y[ans],digits=2),col = "blue", pos = 3, offset = 0.5, font=2, xpd = TRUE)
        if(annotate==2) text(x[ans],y[ans],round(x[ans],digits=2),col = "blue", pos = 1, offset = 1.2, font=2, xpd = TRUE)
        if(annotate==2) text(x[ans],y[ans],round(y[ans],digits=2),col = "blue", pos = 1, offset = 0.5, font=2, xpd = TRUE)
        if(verbose && setIn==1) cat(ans," ", x[ans]," ",y[ans],"\n")
        if(verbose && setIn==2) cat(ans," ",loc[ans]," ",x[ans]," ",y[ans],"\n")
    }
    res
}
   
    if(verbose) 
     {
       if(setIn==1) cat("\nSELECTED DATA POINTS: (index, x, y)\n")
       if(setIn==2) cat("\nSELECTED DATA POINTS: (index, location, x, y)\n")
     }  
    if(setIn==1) pts <- identifyPch(xx,yy,cex=ptsize)
    if(setIn==2) pts <- identifyPch(xx,yy,loc,cex=ptsize)

    if(output==1) 
     { 
       out <- data.frame(cbind(xx[pts],yy[pts]))
       colnames(out) <- c("x","y")
       if(!iso) return(out)
       if(iso && setIn==1) 
        {
          if(verbose) cat("\n * Isolating data between minimum and maximum selected x-values\n")
          out2=iso(cbind(xx,yy),xmin=min(out[,1]),xmax=max(out[,1]),genplot=F,verbose=F)
          return(out2)
        } 
     }  
        
    if(output==2) 
     {
       if(setIn==1) 
        {
          out <- data.frame(cbind(pts,xx[pts],yy[pts]))
          colnames(out) <- c("index","x","y")
        }          
       if(setIn==2)  
        {
          out <- data.frame(cbind(pts,loc[pts],xx[pts],yy[pts]))
          colnames(out) <- c("index","location","x","y")
        }  
       return(out)       
     }  

### END function idPts
}