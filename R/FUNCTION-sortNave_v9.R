### This function is a component of astrochron: An R Package for Astrochronology
### Copyright (C) 2021 Stephen R. Meyers
###
###########################################################################
### sortNave function - (SRM: November 23, 2012, Dec. 2, 2012; May 17, 2013; 
###                           May 20, 2013; June 5, 2013; June 28, 2013;
###                           July 27, 2013; Nov. 27, 2013; June 25, 2014;
###                           March 20, 2017; August 13, 2018; January 14, 2021)
###
### This script will sort data, and calls Fortran routine to
### average duplicates
###########################################################################

sortNave <- function (dat,sortDecr=F,ave=T,xmin=NULL,xmax=NULL,genplot=1,verbose=T)
{

   if(verbose) cat("\n----- PREPARING STRATIGRAPHIC SERIES -----\n")

   dat <- data.frame(dat)
   npts <- length(dat[,1]) 
   if(verbose) cat("\n * Number of stratigraphic samples (rows)=", npts,"\n")
   if(verbose) cat(" * Sorting data, will remove empty entries (from either column).\n")
   
### Remove NAS entries in second column (includes those listed as 'NA')
   dat=subset(dat,!is.na(dat[,2]))  
### sort into increasing or decreasing depth/height/time. Will remove NAS values (in depth column)
   dat <- dat[order(dat[,1],na.last=NA,decreasing=sortDecr),]

   npts <- length(dat[,1])
   if(verbose) cat(" * Number of samples post-sorting=", npts,"\n")

   if(is.null(xmin)) xmin=min(dat[,1])
   if(is.null(xmax)) xmax=max(dat[,1])
   xlim=c(xmin,xmax)
   
### function dup: average duplicates/triplicates/etc.
dup <- function (ipts,x,y)
 {
    F_dat = .Fortran('dupmean_r',
    
    ipts=as.integer(ipts),x=as.double(x),y=as.double(y),
    
    npts=integer(1),xout=double(ipts),yout=double(ipts)
    )

# return the results
    return(cbind(F_dat$xout[1:F_dat$npts],F_dat$yout[1:F_dat$npts]))
  }
### END function dup

### if there are duplicate values, average them (call to dupmean.f)
   t1<-dat[1:(npts-1),1]
   t2<-dat[2:(npts),1]
   dt=t2-t1
   mindt=min(dt)
# save a copy of dat before averaging   
   dat1=dat 
   if(mindt < 1.11022302E-13) 
     {
       cat(" * Duplicates found\n")
       if(ave)
         {
           cat(" * Duplicates values will be averaged.\n")
### call to Fortran routine for quick duplicate averaging
            dat2 <- dup(npts,dat[,1],dat[,2])
            dat2 <- data.frame(dat2)
            npts <- length(dat2[,1]) 
            cat(" * Number of samples post-averaging=",npts,"\n")
            colnames(dat2)[1] <- colnames(dat[1])
            colnames(dat2)[2] <- colnames(dat[2])
            dat <- dat2
          }
       if(!ave) {cat(" * Duplicates found, but will not be averaged.\n")}
     }

### now evaluate sampling statistics
     t1<-dat[1:(npts-1),1]
     t2<-dat[2:(npts),1]
     dt=t2-t1
     dtMin=min(dt) 
     dtMax=max(dt)
     dtMean=mean(dt)     
     dtMedian=median(dt)

     if(verbose)
      {
        cat("\n * Mean sampling interval=", dtMean,"\n")
        cat(" * Median sampling interval=",dtMedian,"\n")
        cat(" * Maximum sampling interval=",dtMax,"\n")
        cat(" * Minimum sampling interval=", dtMin,"\n")
      }
      
### plot data series. Note, cex is the factor by which to increase or decrease default symbol size
     if(genplot==1)
      {
        par(mfrow=c(2,2))
        plot(dat,type="l",col="gray",xlim=xlim,xlab=colnames(dat[1]),ylab=colnames(dat[2]),main="Stratigraphic Series")
        mtext("(red=original data, black=post-averaging)",cex=0.8)
        points(dat1,cex=0.4,col="red")
        points(dat,cex=0.3)
### plot the denisty and the histogram together
        hist(dat[,2],freq=F,xlab=colnames(dat[2]),main="Distribution values"); lines(density(dat[,2], bw="nrd0"),col="red"); grid() 
### boxplot
        boxplot(dat[,2],ylab=colnames(dat[2]))
### Normal probabilty plot (Normal Q-Q Plot)
        qqnorm(dat[,2]); qqline(dat[,2], col="red"); grid()
      }
    
### plot data series. Note, cex is the factor by which to increase or decrease default symbol size
     if(genplot==2)
      {
        par(mfrow=c(1,1))
        plot(dat,xlim=xlim,ylim=c(min(dat1[,2]),max(dat1[,2])),type="l",col="gray",xlab=colnames(dat[1]),ylab=colnames(dat[2]),main="Stratigraphic Series")
        mtext("(red=original data, black=post-averaging)")
        points(dat1,cex=0.4,col="red")
        points(dat,cex=0.3)
      }
        
     return(data.frame(dat))

### END function sortNave
}
