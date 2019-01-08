### This function is a component of astrochron: An R Package for Astrochronology
### Copyright (C) 2018 Stephen R. Meyers
###
###########################################################################
### imbrie: Imbrie and Imbrie (1980) ice model. Model follows convention
###         used in Analyseries.
###          (SRM: April 14-16, 2015; June 19, 2018; December 5-6, 2018; 
###                January 7, 2019)
###
###########################################################################

imbrie <- function (insolation=NULL,Tm=17,b=0.6,times=NULL,initial=0,burnin=100,standardize=T,output=T,genplot=1,verbose=T)
{

  if(verbose) cat("\n----- GENERATING IMBRIE and IMBRIE (1980) ICE MODEL -----\n")
      
  if(is.null(insolation))
   {
     if(verbose) cat(" * Using default insolation series: 65 deg. North, 21 June\n" )
     insolation=getLaskar("insolation",verbose=F)
     insolation=iso(insolation,xmin=0,xmax=1100,genplot=F,verbose=F)
    }  

  if(!is.null(insolation)) insolation=data.frame(insolation)
  npts=length(insolation[,1])

  if(is.null(times)) 
   {
     if(length(Tm) != length(b))
      {
        cat("\n**** ERROR: vectors Tm and b do not have matching entries\n")
        stop("    TERMINATING NOW!")
      }
   }
   
   if(!is.null(times))
    {
      if(length(Tm) != length(times) || length(b) != length(times))
       {
        cat("\n**** ERROR: vectors Tm, b and times do not have matching entries\n")
        stop("    TERMINATING NOW!")
      }
      if(insolation[npts,1] > max(times))
        {
        cat("\n**** ERROR: maximum time for model (in vector times) must be >= maximum time in insolation\n")
        stop("    TERMINATING NOW!")
      }
   }     


# flip so we are going from past to present
  insolation <- insolation[order(insolation[1], na.last = NA, decreasing = T),]
  insolation[1] <- insolation[1]*-1  

# center the insolation to mean of zero
  insolation[2] <- insolation[2]-colMeans(insolation[2])

  y<-double(npts)
  dydt<-double(npts)
  y[1]<-initial
# follow convention of Analyseries  
  x<-insolation[,2]*-1

  if(is.null(times))
   {
# caculate first dydt
# follow convention of Analyseries.
# if present ice sheet index is less than or equal to -insolation, grow  
     if(x[1]>y[1]) dydt[1] <- ( (1-b) / Tm ) * ( x[1]-y[1] )
# otherwise, decay
     if(x[1]<=y[1]) dydt[1] <- ( (1+b) / Tm ) * ( x[1]-y[1] )
# now calculate for remainder of record
     for(i in 2:npts)
      {
        y[i]<-y[i-1] + dydt[i-1]
        if(x[i]>y[i]) dydt[i] <- ( (1-b) / Tm ) * ( x[i]-y[i] )
        if(x[i]<=y[i]) dydt[i] <- ( (1+b) / Tm ) * ( x[i]-y[i] )
      }
   }

  if(!is.null(times))
   {
     ice=data.frame(cbind(times,Tm,b))    
# also flip sequence of ice model parameters to match insolation direction
     ice <- ice[order(ice[1], na.last = NA, decreasing = T),]
     ice[1] <- ice[1]*-1
# find first value to use
     getTime=max(which( ice[1] <= insolation[1,1] ))
     if(x[1]>y[1]) dydt[1] <- ( (1-ice[getTime,3]) / ice[getTime,2] ) * ( x[1]-y[1] )
# otherwise, decay
     if(x[1]<=y[1]) dydt[1] <- ( (1+ice[getTime,3]) / ice[getTime,2] ) * ( x[1]-y[1] )
# now calculate for remainder of record
     for(i in 2:npts)
      {
        getTime=max(which( ice[1] <= insolation[i,1] ))
        y[i]<-y[i-1] + dydt[i-1]
        if(x[i]>y[i]) dydt[i] <- ( (1-ice[getTime,3]) / ice[getTime,2] ) * ( x[i]-y[i] )
        if(x[i]<=y[i]) dydt[i] <- ( (1+ice[getTime,3]) / ice[getTime,2] ) * ( x[i]-y[i] )
      }
    }

# remove burn in points, and make a data frame
  out<-data.frame(cbind(insolation[burnin:npts,1],y[burnin:npts]))

# flip to original direction
  out <- out[order(out[1], na.last = NA, decreasing = T),]
  out[1] <- out[1]*-1  

  if(standardize) 
   {
     out[2] <- out[2]-min(out[2])
     out[2] <- out[2]/max(out[2])
   }  

if(genplot==1 || genplot==2)
  {
# remove burn in points, and make a data frame
  insolation<-data.frame(cbind(insolation[burnin:npts,1],insolation[burnin:npts,2]))

# flip to original direction
  insolation <- insolation[order(insolation[1], na.last = NA, decreasing = T),]
  insolation[1] <- insolation[1]*-1  
  
### plots
   if(genplot==1) 
    {
     par(mfrow=c(2,1))
     plot(insolation,cex=0.5,xlab="Time (ka BP)",ylab="Insolation",main="Insolation model",type="l",col="red")
#    plot(out,cex=0.5,xlab="Time (ka BP)",ylab="Ice volume",main="Ice volume model",type="l",ylim=c(max(out[,2]),min(out[,2])))
     plot(out,cex=0.5,xlab="Time (ka BP)",ylab="Ice volume",main="Ice volume model",type="l",ylim=c(min(out[,2]),max(out[,2])))
     polygon(c(min(out[,1]),out[,1],max(out[,1])),c(min(out[,2]),out[,2],min(out[,2])),col="cornflowerblue",border=NA)
   }

   if(genplot==2) 
    {
      dev.new(height=6,width=10)
      layout(matrix(c(1,2,3,3), 2, 2, byrow = FALSE))
      pts2=length(out[,2])
      xlim=c(min(insolation[1:pts2,2]),max(insolation[1:pts2,2]))
      ylim=c(min(out[1:pts2,2]),max(out[1:pts2,2]))
      for(i in pts2:1)
       {
         plot(insolation[i:pts2,1],insolation[i:pts2,2],cex=0.5,xlab="Time (ka BP)",ylab="Insolation",main="Insolation model",type="l",ylim=xlim,col="red")
         plot(out[i:pts2,1],out[i:pts2,2],cex=0.5,xlab="Time (ka BP)",ylab="Ice volume",main="Ice volume model",type="l",ylim=ylim,col="cornflowerblue")
         plot(insolation[i:pts2,2],out[i:pts2,2],type="l",xlab="Insolation",ylab="Ice volume",col="darkgreen",lty=2,xlim=xlim,ylim=ylim)
         points(insolation[i,2],out[i,2],col="darkgreen",pch=16)
         Sys.sleep(0.15)
       }
    }    
  }
  
  if(output) return(out)

### END function imbrie
}
