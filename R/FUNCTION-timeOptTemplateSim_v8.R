### This function is a component of astrochron: An R Package for Astrochronology
### Copyright (C) 2021 Stephen R. Meyers
###
###########################################################################
### function timeOptTemplateSim - (SRM: May 28, 2012; Oct. 14, 2014; Oct. 17, 2014;
###                             Oct. 19, 2014; Jan. 13, 2015; March 9, 2015
###                             June 8, 2015; Sept. 30, 2015; 
###                             October 20-21, 2015; November 19, 2015;
###                             December 17, 2015; February 7, 2016; 
###                             February 25, 2016; October 18-26, 2016;
###                             December 13-21, 2017; January 1, 2018; 
###                             November 24, 2018; January 14, 2021)
###########################################################################

timeOptTemplateSim <- function (dat,template=NULL,corVal=NULL,numsim=2000,rho=NULL,sedmin=0.5,sedmax=5,difmin=NULL,difmax=NULL,fac=NULL,numsed=50,linLog=1,limit=T,fit=1,fitModPwr=T,iopt=3,flow=NULL,fhigh=NULL,roll=NULL,targetE=NULL,targetP=NULL,cormethod=1,detrend=T,detrendTemplate=F,flipTemplate=F,ncores=1,output=0,genplot=T,check=T,verbose=T)
{

if(verbose) cat("\n----- TimeOptTemplate Monte Carlo Simulation -----\n")

cormethod=1

# prepare data array
   dat = data.frame(dat)      
   npts <- length(dat[,1]) 
   dx <- dat[2,1]-dat[1,1]

if(check)
 {
# error checking 
   if(dx<0)
     { 
       if (verbose) cat("\n * Sorting data into increasing height/depth/time, removing empty entries\n")
       dat <- dat[order(dat[,1], na.last = NA, decreasing = F), ]
       dx <- dat[2,1]-dat[1,1]
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
   cat(" * Stratigraphic series length (meters):",(npts-1)*dx,"\n")
   cat(" * Sampling interval (meters):",dx,"\n\n")
 }

# if the correlation value from the data analysis wasn't entered, calculate it now.
   if(is.null(corVal))
    {
# parallelization will be activated in timeOptTemplate
     res=timeOptTemplate(dat,template=template,sedmin=sedmin,sedmax=sedmax,difmin=difmin,difmax=difmax,fac=fac,numsed=numsed,linLog=linLog,limit=limit,fit=fit,fitModPwr=fitModPwr,iopt=iopt,flow=flow,fhigh=fhigh,roll=roll,targetE=targetE,targetP=targetP,cormethod=cormethod,detrend=detrend,detrendTemplate=detrendTemplate,flipTemplate=flipTemplate,ncores=ncores,check=T,output=1,genplot=1,verbose=1)
     datCorPwr = max(res[,2])
    }
   if(!is.null(corVal)) datCorPwr=corVal
   
if(verbose)
 {  
  cat(" * (Envelope r^2) x (Spectral Power r^2) =", datCorPwr,"\n")
 }

#######################################################################################
# Monte Carlo Simulation

# calculate the ar1 coeff for the simulations. be sure to prepare data set as
#  was done in timeOptTemplate
# detrend
if (detrend) 
  {
    lm.1 <- lm(dat[,2] ~ dat[,1])
    dat[2] <- dat[2] - (lm.1$coeff[2]*dat[1] + lm.1$coeff[1])
    if(verbose) cat(" * Linear trend subtracted. m=",lm.1$coeff[2],"b=",lm.1$coeff[1],"\n")
  }

# standardize data series
   dat[2]=dat[2]-colMeans(dat[2])
   dat[2]=dat[2]/sapply(dat[2],sd)

### what is the estimated AR1 coefficient?
   if(is.null(rho))
    {
      lag0 <- dat[1:(npts-1),2]
      lag1 <- dat[2:npts,2]
      rho <- cor(lag0,lag1)
      if(verbose) cat(" * Raw AR1 =",rho,"\n")
    }  

if(verbose) 
  {
    cat("\n * PLEASE WAIT: Performing",numsim,"Simulations \n")
    cat("\n0%       25%       50%       75%       100%\n")
# create a progress bar
    progress = utils::txtProgressBar(min = 0, max = numsim, style = 1, width=43)
  }


#  create output array, dimension appropriately
#  simres will contain envelope*power
simres <- double(numsim)

# begin simulation loop
isim=0
for (isim in 1:numsim) 
  {
    if(verbose) utils::setTxtProgressBar(progress, isim)
# generate AR1 noise
    sim = ar1(npts, dx, mean=0, sdev=1, rho=rho, genplot=F, verbose=F)
    sim[1]=dat[1]
# recenter and standardize
    sim[2]=sim[2]-colMeans(sim[2])
    sim[2]=sim[2]/sapply(sim[2],sd)               
    simres[isim]=max(timeOptTemplate(sim,template=template,sedmin=sedmin,sedmax=sedmax,difmin=difmin,difmax=difmax,fac=fac,numsed=numsed,linLog=linLog,limit=limit,fit=fit,iopt=iopt,flow=flow,fhigh=fhigh,roll=roll,targetE=targetE,targetP=targetP,cormethod=cormethod,detrend=detrend,detrendTemplate=detrendTemplate,flipTemplate=flipTemplate,ncores=ncores,check=F,output=1,genplot=0,verbose=0)[2])
    if(verbose) cat("Simulation",isim,"r2=",simres[isim],"\n")
# end simulation loop
}

if(verbose) close(progress) 

# now sort results, determine how many have values > your result
# envelope * spectral power
    numgt = sum(simres>datCorPwr)
    pvalCorPwr=numgt/numsim 
    if(pvalCorPwr < (10/numsim) && (10/numsim) <=1 ) pvalCorPwr= 10/numsim
    if(pvalCorPwr >= (10/numsim) && (10/numsim) <=1 ) pvalCorPwr=pvalCorPwr   
    if((10/numsim) > 1 ) pvalCorPwr=1
     
    if(verbose) cat("\n * (Envelope r^2) * (Spectral Power r^2) p-value =",pvalCorPwr, "\n")

     
    if(genplot)
     { 
      dev.new(title = paste("TimeOpt Monte Carlo Results"), height = 5, width = 6)
      par(mfrow=c(1,1))
      plot(density(simres), col="black",xlim=c(0,1),type="l",xlab=expression(paste({"r"^2}["opt"])),main=expression(bold(paste({"r"^2}["opt"]," Monte Carlo Results"))),cex.lab=1.1,lwd=2)
      polygon(density(simres),col="red",border=NA)
#      grid()
      abline(v=datCorPwr,col="blue",lwd=2,lty=3)
      mtext(round(datCorPwr,digits=5),side=3,line=0,at=datCorPwr,cex=1,font=4,col="blue")
     }
     
# output = (0) nothing, (1) envelope*spectral power r^2 p-value, (2) output simulation r^2 results
     if(output == 1) return(data.frame(pvalCorPwr))
     if(output == 2) return(data.frame(simres))
     
### END function timeOptTemplateSim
}
