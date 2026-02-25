### This function is a component of astrochron: An R Package for Astrochronology
### Copyright (C) 2026 Stephen R. Meyers
###
##################################################################################
### function timeOptBMCMC - (SRM: Feb. 23, 2026)
##################################################################################

timeOptBMCMC <- function(dat,age=NULL,env=TRUE,mvopt=0,sedmin=0.5,sedmax=5,ftol=0.005,
                         roll=10^7,nsamples=50000,test=FALSE,detrend=TRUE,
                         savefile=FALSE,genplot=TRUE,check=TRUE,verbose=TRUE)
 {
  if(verbose) 
   {
    cat("\n----- TimeOptBMCMC: TimeOpt Bayesian Astrochronologic Evaluation -----\n")
    cat("             Using the Method of Malinverno and Meyers (2024) \n\n")
   }
  if(test) cat("**** WARNING: testing mode activated. All MCMC samples accepted!\n\n")

#######################################################################################
#### (1) data preparation
#######################################################################################
# prepare data array
  dat = data.frame(dat)
# number of data points   
  Ndv = length(dat[,1]) 
  dx = dat[2,1]-dat[1,1]

# error checking 
  if(check)
   {
    if(dx<0)
     { 
      if (verbose) cat("\n * Sorting data into increasing height/depth/time, removing empty entries\n")
      dat = dat[order(dat[1], na.last = NA, decreasing = F), ]
      dx = dat[2,1]-dat[1,1]
      Ndv = length(dat[,1])
     }
    dtest = dat[2:Ndv,1]-dat[1:(Ndv-1),1] 
    epsm = 1e-9
    if( (max(dtest)-min(dtest)) > epsm ) cat("\n**** ERROR: sampling interval is not uniform. TERMINATING NOW!\n")
    if(is.null(age)) stop("\n**** ERROR: age has not been specified. TERMINATING NOW!\n")
    if(age <= 0) stop("\n**** ERROR: age must be > 0. TERMINATING NOW!\n")
    if(nsamples < 50000) cat("\n**** WARNING: nsamples should typically be at least 50000\n\n\n")
   }
  if (verbose) 
   {
    cat(" * Number of data points in stratigraphic series:",Ndv,"\n")
    cat(" * Stratigraphic series length (meters):",(Ndv-1)*dx,"\n")
    cat(" * Sampling interval (meters):",dx,"\n")
   }
  if (detrend) 
   {
    lm.1 <- lm(dat[,2] ~ dat[,1])
    dat[2] = dat[2] - (lm.1$coeff[2]*dat[1] + lm.1$coeff[1])
    if(verbose) cat(" * Linear trend subtracted. m=",lm.1$coeff[2],"b=",lm.1$coeff[1],"\n")
   }
# standardize data series
  dat[2] = dat[2]-colMeans(dat[2])
  dat[2] = dat[2]/sapply(dat[2],sd)
  
# convert sedmin and sedmax from cm/kyr to m/kyr
  sedmin = sedmin/100
  sedmax = sedmax/100 

  fconv = 1000/(360*60*60)             # conversion from arcsec/yr to cycles/kyr
  
#######################################################################################
#### (2) if check=T, ensure resolution requirements are satisfied
#######################################################################################
  if(check)
   {
    prior <- .getpriorpdfastropara(age/1000,mvopt)             # parameters for priors, excluding u
# check minimum and maximum sedimentation rates. sedmin is in meters/kyr, dx is in meters.
    NyqFreq=sedmin/(2*dx)
# freqHigh identifies the frequency of the shortest period in the target, given the 
#   stated uncertainties in g and k. we add ftol for the maximum possible upper half power 
#   point in the taner filtering.  we will take 2 standard deviations of mean as limit.
    eop <- .compeopfreq(c(prior$snmean+2*prior$sigma,0),mvopt)
    pvP=eop$pv*fconv
    freqHigh=max(pvP)+ftol
    if(freqHigh>NyqFreq)
     {
      sedmin = 2*dx*freqHigh
      if(verbose) cat("\n**** WARNING: minimum sedimentation rate is too low.\n")
      if(verbose) cat("            -> sedmin reset to",100*sedmin,"cm/ka\n")
     }
# check maximum sedimentation rate. sedmax is in meters/kyr. dx is in meters.
    RayFreq = sedmax/(Ndv*dx)
# freqLow identifies the frequency of the longest period in the target, given the 
#   stated uncertainties in g. we will not worry about ftol here. we will take 2 standard 
#   deviations of mean as nominal limit.
    eop <- .compeopfreq(c(prior$snmean-2*prior$sigma,0),mvopt)
    evP=eop$ev*fconv
    freqLow=min(evP)
    if(RayFreq>freqLow)
     {
      sedmax = Ndv*dx*freqLow
      if(verbose) cat("\n**** WARNING: maximum sedimentation rate is too high.\n")
      if(verbose) cat("            -> sedmax reset to",100*sedmax,"cm/ka\n")
     }
    if(sedmin>sedmax) stop("\n**** ERROR: sedmin > sedmax. TERMINATING NOW!")
   }

#######################################################################################
#### (3) set parameters in pdfpara, and additional variables/constants
#######################################################################################
# pdfpara 'structure' passed to mcmcadapt and other functions. 
# pdfpara will accompany:
#  dat (data frame) 
#    - depth 'vector' (zv of MATLAB)
#    - data 'vector' (dv of MATLAB)
#  fc_dat (data frame)
#    - 'vector' of spatial frequencies
#    - 'vector' of real fourier coeff of dat (pdfpara.ftdv of MATLAB)
#    - 'vector' of imaginary fourier coeff of dat (pdfpara.ftdv of MATLAB)
#  mv (model vector)
#    - The first 5 elements of mv are always g_1 to g_5
#    - The last 2 elements of mv are always k and u, respectively
#    - The first 10 columns of Gdv will contain sine and cosine terms for 5
#      eccentricity freqs. that are differences g_i - g_j
#  prior (list)
#    - parameters defining the priors, derived from getpriorpdfastropara 
#       (mu,sigma,alpha,a1,mu1,sigma1,snmean)
  pdfpara = double(14)
  names(pdfpara) = c('tGa','prioronly','P','roll','deltafp','likfonly','nolikfdenv',
                     'mvopt','iu','umin','umax','ik','taudenv','ar2')
  pdfpara['tGa'] = age/1000            # age in billions of years
  if(test) pdfpara['prioronly'] = 1    # exclude likelihood, only evaluate prior for testing
  if(!test) pdfpara['prioronly'] = 0   # conduct complete analysis
  pdfpara['P'] = 2                     # order of AR(P) process fitted to the residuals
  pdfpara['roll'] = roll               # roll-off rate in Taner filter (dimensionless)
  pdfpara['deltafp'] = ftol            # default delta frequency for Taner filter boundaries
  pdfpara['likfonly'] = 0              # (default 0) if 1 use the likelihood as the target pdf
  if(env) pdfpara['nolikfdenv'] = 0    # (default 0) if zero use the envelope likelihood
  if(!env) pdfpara['nolikfdenv'] = 1   # omit the envelope likelihood
  pdfpara['mvopt'] = mvopt             # variable that defines the astro. parameters in mv
  pdfpara['iu'] = .getmvindex('u',mvopt) # index of u (sedimentation rate) in mv
  pdfpara['umin'] = sedmin             # minimum sedimentation rate to evaluate (m/kyr)
  pdfpara['umax'] = sedmax             # maximum sedimentation rate to evaluate (m/kyr)
  pdfpara['ik'] = .getmvindex('k',mvopt) # index of k (precession frequency) in mv
  pdfpara['ar2'] = 4                   # option for ar.burg model fitting
  
# taudenv will be added later as pdfpara[14]
# note: the vector 'zv' of the MATLAB pdfpara is carried along in dat[1]
# note: the vector 'dv' of MATLAB is in dat[2]

# calculate fourier coefficients for the stratigraphic series. the dataframe fc_dat
#  contains three columns: frequency, real fourier coeff, imaginary fourier coeff
  fc_dat <- periodogram(dat,padfac=2,demean=T,detrend=F,output=2,nrm=0,genplot=F,check=F,verbose=F)

# set additional variables/constants
  niterinbatch = 50                    # adaptive MCNC: number of samples per batch
  if(!check) prior <- .getpriorpdfastropara(pdfpara['tGa'],pdfpara['mvopt'])  # parameters for priors, excluding u

#######################################################################################
#### (4) conduct Bayesian inversion
#######################################################################################
# set starting value of parameter vector to prior means for g_i, s_i
  priormeanastro = prior$mu
  priorsdevastro = prior$sigma
  nmv = length(priormeanastro)+1      # g's, s's, k, and u
  mvstart = double(nmv)
  mvstart[1:(nmv-1)] = priormeanastro

# for k and u, set the starting value parameter vector to a slightly perturbed value from
# the prior pdf, so each sampling chain starts from a different point
  ik = pdfpara['ik']                   # index to k in mv
  iu = pdfpara['iu']                   # index to u in mv
  kstart = priormeanastro[ik]+0.1*priorsdevastro[ik]*rnorm(1)  # starting value of k
  mvstart[ik] = kstart
  umid = (pdfpara['umin']+pdfpara['umax'])/2
  urange = pdfpara['umax']-pdfpara['umin']
  ustart <- runif(1,min=umid-0.1*urange,max=umid+0.1*urange)   # starting value of sed.rate
  mvstart[iu] = ustart

# set initial standard deviations for candidate generation
  sigmamvcand = double(nmv)
  sigmamvcand[1:(nmv-1)] = priorsdevastro[1:(nmv-1)]
  sigmamvcand[ik] = min(0.5,priorsdevastro[ik])                # prior st. dev of k can be large
  sigmamvcand[iu] = 0.01*ustart

# compute period of shortest eccentricity cycle and taudenv
  if (pdfpara['prioronly'] == 1) pdfpara['taudenv']= NA
  if (pdfpara['prioronly'] == 0)
   {
    ev <- .compeopfreq(mvstart,pdfpara['mvopt'])$ev
    periode=1/(fconv*max(ev))          # period (kyr) of shortest eccentricity cycle
# tau for envelope residuals at u=umid
    pdfpara['taudenv']=umid*periode/(4*dx)
    cat('\n * taudenv:', pdfpara['taudenv'],"\n")
   }

# prepare output files if desired    
  if(savefile) 
   {
    mvlabels <- .getmvlabels(pdfpara['mvopt'])     
    write.table(file="mcmc-mv.csv",matrix(data=c('iter','logpdf',mvlabels),nrow=1,ncol=length(mvlabels)+2), 
                sep = ",", row.names = FALSE, col.names= FALSE, append=FALSE) 
    write.table(file="mcmc-pacc.csv", matrix(data=c('iter',mvlabels),nrow=1,ncol=length(mvlabels)+1), 
                sep = ",", row.names = FALSE, col.names= FALSE, append=FALSE)   
    write.table(file="mcmc-sigmamvcand.csv", matrix(data=c('iter',mvlabels),nrow=1,ncol=length(mvlabels)+1), 
                sep = ",", row.names = FALSE, col.names= FALSE, append=FALSE)      
   }

# ===== run mcmcadapt() =====   
  mcmcres <- .mcmcadapt(mvstart,dat,fc_dat,prior,pdfpara,sigmamvcand,niter=nsamples,niterinbatch,savefile)
  
# ===== postprocess and plot =====
  if(!test)
   {
# burn in detection, using median value from second half of sample pdf
    logpdfburnin=median(mcmcres$mv[(nsamples/2):nsamples,1])
# identify the first value to above it, use that as the start of the burn-in  
    nburnin=which(mcmcres$mv[,1]>logpdfburnin)[1]
    if(verbose) cat(" * MCMC chain burn-in: discard",nburnin-1, "samples.\n")
   }else{
    nburnin=1
   }

# diagnostic plots and summary statistics
  if(genplot) timeOptBMCMCplot(dat,mcmcres,pdfpara,nburnin,fc_dat)

  return(invisible(list(mcmcres=mcmcres,pdfpara=pdfpara,dat=dat)))
# end function timeOptBMCMC
}  
