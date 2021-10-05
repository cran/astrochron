### This function is a component of astro: An R Package for Astrochronology
### Copyright (C) 2021 Stephen R. Meyers
###
###########################################################################
### eTimeOpt - (SRM: February 26-28, 2016; October 5, 2016; 
###                  October 13, 2017; December 7, 2017; December 18, 2017
###                  December 4, 2021; January 21, 2021; August 30, 2021)
###
### evolutive version of timeOpt
###
### this version is a wrapper that loops timeOpt for evolutive analysis
###########################################################################

eTimeOpt <- function (dat,win=dt*100,step=dt*10,sedmin=0.5,sedmax=5,numsed=100,linLog=1,limit=T,fit=1,fitModPwr=T,flow=NULL,fhigh=NULL,roll=NULL,targetE=NULL,targetP=NULL,detrend=T,ydir=1,palette=6,ncolors=100,output=1,genplot=T,check=T,verbose=1)
{

if(verbose==0 || verbose==1 || verbose==2) timeOptVerbose=F
if(verbose==3) timeOptVerbose=T

# plotting library
#library(fields)

if(verbose==1 || verbose==2 || verbose==3) cat("\n----- PERFORMING EVOLUTIVE TimeOpt ANALYSIS -----\n")

dat <- data.frame(dat)
npts <- length(dat[,1]) 
dt <- dat[2,1]-dat[1,1]
 
if(check)
 {
# error checking 
   if(dt<0)
     { 
       if (verbose==1 || verbose==2 || verbose==3) cat("\n * Sorting data into increasing height/depth/time, removing empty entries\n")
       dat <- dat[order(dat[,1], na.last = NA, decreasing = F), ]
       dt <- dat[2,1]-dat[1,1]
       npts <- length(dat[,1])
     }
   dtest <- dat[2:npts,1]-dat[1:(npts-1),1] 
   epsm=1e-9
   if( (max(dtest)-min(dtest)) > epsm ) 
     {
       cat("\n**** ERROR: dat sampling interval is not uniform.\n")
       stop("**** TERMINATING NOW!")
     }
}

   if( palette != 1 && palette != 2 && palette != 3 && palette != 4 && palette != 5 && palette != 6) 
     {
       cat("\n**** WARNING: palette option not valid. Will use palette = 6.\n")
       palette = 6
     }

if (verbose ==1 || verbose ==2 || verbose ==3) 
 {
   cat(" * Number of data points in stratigraphic series:",npts,"\n")
   cat(" * Stratigraphic series length (space or time):",(npts-1)*dt,"\n")
   cat(" * Sampling interval (space or time):",dt,"\n")
 }

# number of points per window
winpts=as.integer( floor(win/dt) + 1 )
# number of points to increment window
wininc=as.integer( floor(step/dt) )
# number of spectra
nspec= floor( (npts-winpts) / wininc) + 1 

# error checking 
if(winpts >= npts) stop(" * Specified duration longer than permitted for evolutive analysis.\n") 

if (verbose==1 || verbose==2 || verbose==3)
 {
   cat("\n * Number of data points per window:",winpts,"\n")
   cat(" * Moving window size (space or time):",(winpts-1)*dt,"\n")
   cat(" * Window step points:",wininc,"\n")
   cat(" * Window step (space or time):",wininc*dt,"\n")
   cat(" * Number of windows:",nspec,"\n")
 }
 

# if sedmin equals sedmax, ensure numsed = 1
  if(sedmin == sedmax) 
   {
     if(numsed != 1)
      {
        if(verbose) cat("\n**** WARNING: sedmin = sedmax, so numsed will be reset to 1.\n\n")
        numsed=1
      }  
   }  

# set up matrices for results
sedrates<-double(numsed)
height<-double(nspec)
r2_env<-double(numsed*nspec)
dim(r2_env)<-c(numsed,nspec)
r2_pwr<-double(numsed*nspec)
dim(r2_pwr)<-c(numsed,nspec)
r2_opt<-double(numsed*nspec)
dim(r2_opt)<-c(numsed,nspec)

if(verbose==1)
 {
   cat("\n * PLEASE WAIT: Performing Optimizations\n")
   cat("\n0%       25%       50%       75%       100%\n")
# create a progress bar
   progress = utils::txtProgressBar(min = 0, max = nspec, style = 1, width=43)
 }

# begin windowing loop
ihold = 1
for (i in 1:nspec)
 {
   if(verbose==1) utils::setTxtProgressBar(progress, i)
   if (verbose ==2 || verbose ==3) 
    {
      cat("\n")
      cat("***********************************************************************************\n")
      cat("  - Processing",i,"of",nspec,"; Points=",ihold,"-",(ihold+winpts-1),"; Mid-Depth=",(dat[ihold,1]+dat[(ihold+winpts-1),1])/2,":", dat[ihold,1],"-",dat[(ihold+winpts-1),1],"\n")
      cat("***********************************************************************************\n")
    }
# prep data
# shuffle in data
   dat2 <- data.frame(cbind( dat[ihold:(ihold+winpts-1),1], dat[ihold:(ihold+winpts-1),2] ))
 
   res <- timeOpt(dat2,sedmin=sedmin,sedmax=sedmax,numsed=numsed,limit=limit,linLog=linLog,fit=fit,fitModPwr=fitModPwr,flow=flow,fhigh=fhigh,roll=roll,targetE=targetE,targetP=targetP,detrend=detrend,output=1,genplot=F,check=F,verbose=timeOptVerbose)
   
### save frequency and power to freq and pwr
   if(i == 1) sedrates <- res[,1]
   r2_env[,i] <- res[,2]
   r2_pwr[,i] <- res[,3]
   r2_opt[,i] <- res[,4]
   height[i]= dat[1,1] + ( (i-1)*dt*wininc ) + ( (winpts-1)*dt/2 )

   ihold = ihold + wininc   
# end windowing loop
}

 if(verbose==1) close(progress)

# now add plots
 if(genplot && numsed>1)
  {

# set color palette
#  rainbow colors
    if(palette == 1) colPalette = tim.colors(ncolors)
#  grayscale
    if(palette == 2) colPalette = gray.colors(n=ncolors,start=1,end=0,gamma=1.75)
#  dark blue scale (from larry.colors)
    if(palette == 3) colPalette = colorRampPalette(c("white","royalblue"))(ncolors)
#  red scale
    if(palette == 4) colPalette = colorRampPalette(c("white","red2"))(ncolors)
#  blue to red plot
    if(palette == 5) colPalette = append(colorRampPalette(c("royalblue","white"))(ncolors/2),colorRampPalette(c("white","red2"))(ncolors/2))
# viridis colormap
    if(palette == 6) colPalette = viridis(ncolors, alpha = 1, begin = 0, end = 1, direction = 1, option = "D")

    dev.new(title=paste("eTimeOpt results"),height=5.3,width=10)
    par(mfrow=c(1,3))
    par(mar=c(4.1, 4.1, 4.1, 5.1))

    if (ydir == 1) 
     {
       image.plot(sedrates,height,r2_env,col = colPalette,useRaster=F,xlab="Sedimentation Rate (cm/ka)",ylab="Height (m)",main=expression(paste("Envelope (r"^"2",""["envelope"],")")),cex.lab=1.3,cex.main=1.4)
       image.plot(sedrates,height,r2_pwr,col = colPalette,useRaster=F,xlab="Sedimentation Rate (cm/ka)",ylab="Height (m)",main=expression(paste("Power (r"^"2",""["power"],")")),cex.lab=1.3,cex.main=1.4)
       image.plot(sedrates,height,r2_opt,col = colPalette,useRaster=F,xlab="Sedimentation Rate (cm/ka)",ylab="Height (m)",main=expression(paste("Envelope*Power (r"^"2",""["opt"],")")),cex.lab=1.3,cex.main=1.4)
     }  
       
    if (ydir == -1) 
     {
# in this case, reset ylim range.
       ylimset=c( max(height),min(height) )
       image.plot(sedrates,height,r2_env,ylim=ylimset,col = colPalette,xlab="Sedimentation Rate (cm/ka)",ylab="Depth (m)",main=expression(paste("Envelope (r"^"2",""["envelope"],")")),cex.lab=1.3,cex.main=1.4)
       image.plot(sedrates,height,r2_pwr,ylim=ylimset,col = colPalette,xlab="Sedimentation Rate (cm/ka)",ylab="Depth (m)",main=expression(paste("Power (r"^"2",""["power"],")")),cex.lab=1.3,cex.main=1.4)
       image.plot(sedrates,height,r2_opt,ylim=ylimset,col = colPalette,xlab="Sedimentation Rate (cm/ka)",ylab="Depth (m)",main=expression(paste("Envelope*Power (r"^"2",""["opt"],")")),cex.lab=1.3,cex.main=1.4)
      }
  }

 if(genplot && numsed==1)
  {
    dev.new(height=7,width=7)
    par(mfrow = c(1, 3))
     
    if (ydir == 1) 
     {
       plot(r2_env,height,type="l",main=expression(paste("Envelope (r"^"2",""["envelope"],")")),xlab="",ylab="Location")
       plot(r2_pwr,height,type="l",main=expression(paste("Power (r"^"2",""["power"],")")),xlab="",ylab="Location")
       plot(r2_opt,height,type="l",main=expression(paste("Envelope*Power (r"^"2",""["opt"],")")),xlab="",ylab="Location")
     }   
    if (ydir == -1) 
     {
# in this case, reset ylim range.
       ylimset=c( max(height),min(height) )   
       plot(r2_env,height,type="l",ylim=ylimset,main=expression(paste("Envelope (r"^"2",""["envelope"],")")),xlab="",ylab="Location")
       plot(r2_pwr,height,type="l",ylim=ylimset,main=expression(paste("Power (r"^"2",""["power"],")")),xlab="",ylab="Location")
       plot(r2_opt,height,type="l",ylim=ylimset,main=expression(paste("Envelope*Power (r"^"2",""["opt"],")")),xlab="",ylab="Location")
     }
  }

# OUTPUT RESULTS
# add column titles to identify each record for output
   colnames(r2_env) <- height
   colnames(r2_pwr) <- height
   colnames(r2_opt) <- height

  if(numsed==1)
   {   
    if (output==1) return( data.frame( cbind(height,r2_env[1,],r2_pwr[1,],r2_opt[1,]) ))
    if (output==2) return( data.frame( cbind(height,r2_env[1,]) ))
    if (output==3) return( data.frame( cbind(height,r2_pwr[1,]) ) )
    if (output==4) return( data.frame( cbind(height,r2_opt[1,]) ) )
   }

  if(numsed>1)
   {   
    if (output==1) 
     {
# add frequency column
       r2_env_out <- data.frame( cbind(sedrates,r2_env) )
       r2_pwr_out <- data.frame( cbind(sedrates,r2_pwr) )
       r2_opt_out <- data.frame( cbind(sedrates,r2_opt) )
       return( list(r2_env_out,r2_pwr_out,r2_opt_out) )
     }

    if (output==2) return( data.frame( cbind(sedrates,r2_env) ))
    if (output==3) return( data.frame( cbind(sedrates,r2_pwr) ) )
    if (output==4) return( data.frame( cbind(sedrates,r2_opt) ) )
   }


#### END function eTimeOpt
}
