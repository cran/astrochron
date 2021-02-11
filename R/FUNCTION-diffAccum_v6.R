### This function is a component of astrochron: An R Package for Astrochronology
### Copyright (C) 2021 Stephen R. Meyers
###
###########################################################################
### diffAccum: model differential accumulation
###              (SRM: September 5-7, 2015; Sept. 21, 2015; 
###                    December 17, 2017; June 12, 2018; January 14, 2021)
###                                        
###########################################################################

# could modify for explicit differential unstretching, if depth input instead of time

diffAccum <- function (dat,sedmin=0.01,sedmax=0.02,dir=1,genplot=T,verbose=T)
{

  if (verbose) cat("\n----- CONVERTING TEMPORAL RECORD TO SPACE USING DIFFERENTIAL ACCUMULATION -----\n")
  dat<-data.frame(dat)
  dt <- dat[2,1]-dat[1,1]
  npts <- length(dat[,1]) 

# error checking 
   if(dt<0)
     { 
       if (verbose) cat("\n * Sorting data into increasing time, removing empty entries\n")
       dat <- dat[order(dat[,1], na.last = NA, decreasing = F), ]
       dt <- dat[2,1]-dat[1,1]
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
   cat(" * Number of data points in series:",npts,"\n")
   cat(" * Series length (time):",(npts-1)*dt,"\n")
   cat(" * Sampling interval (time):",dt,"\n")
 }

# convert input variable to range from sedmin to sedmax
  x=dat[,2]
  if(dir == 1) 
   {
     slope=(sedmax-sedmin)/(max(x)-min(x))
     b=sedmax-(slope*max(x))
   }
  if(dir == -1) 
   {
     slope=(sedmax-sedmin)/(min(x)-max(x))
     b=sedmin-(slope*max(x))
   }
   
  sedrate=(slope*x)+b
  if (verbose) cat(" * Mean sedimentation rate=:",mean(sedrate),"m/ka \n")

# model differential accumulation
# ensure time starts at zero
  out<-dat
  out[,1]<-dat[,1]-dat[1,1]
  for (i in 2:npts) 
   {
      out[i,1]=out[(i-1),1]+dt*sedrate[i]
   }

if(genplot)
  {
### plots
    par(mfrow=c(3,1))
    plot(dat,cex=0.5,xlab="Time (ka)",ylab="Value",main="Times Series"); lines(dat)
    plot(dat[,1],sedrate,cex=0.5,xlab="Time (ka)",ylab="Sedimentation Rate",main="Sedimentation Model"); lines(dat[,1],sedrate)
    plot(out,cex=0.5,xlab="Depth (meters)",ylab="Value",main="Modeled Stratigraphic Series"); lines(out)

### plot the denisty and the histogram together
#    hist(dat3[,2],freq=F); lines(density(dat3[,2], bw="nrd"),col="red"); grid()
### boxplot
#    boxplot(dat3[,2])
### Normal probabilty plot (Normal Q-Q Plot)
#    qqnorm(dat3[,2]); qqline(dat3[,2], col="red")
  }
  
  return(out)

### END function diffAccum
}
