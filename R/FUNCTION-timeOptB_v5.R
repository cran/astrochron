### This function is a component of astrochron: An R Package for Astrochronology
### Copyright (C) 2026 Stephen R. Meyers
###
###########################################################################
### function timeOptB - (SRM: Feb. 24, 2026)
###########################################################################
timeOptB <- function (dat,age=NULL,env=TRUE,mvopt=0,sedmin=0.5,sedmax=5,numsed=120,
            kmin=NULL,kmax=NULL,numk=120,ftol=0.005,roll=10^7,detrend=TRUE,
            palette=6,genplot=TRUE,check=TRUE,verbose=TRUE)
 {

  if(verbose) 
   {
    cat("\n----- TimeOptB: TimeOpt Bayesian Astrochronologic Evaluation -----\n")
    cat("          Using the Method of Malinverno and Meyers (2024) \n\n")
   }

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
    if( (max(dtest)-min(dtest)) > epsm ) stop("\n**** ERROR: sampling interval is not uniform. TERMINATING NOW!\n")
    if(is.null(age)) stop("\n**** ERROR: age has not been specified. TERMINATING NOW!\n")
    if(age <= 0) stop("\n**** ERROR: age must be > 0. TERMINATING NOW!\n")
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
  umin = sedmin/100
  umax = sedmax/100 

  fconv = 1000/(360*60*60)             # conversion from arcsec/yr to cycles/kyr
  k0=50.475838                         # present precession freq. (arcsec/yr), Laskar et al. 2004
  tGa=age/1000                         # age in billions of years
  Nu=numsed                            # number of sedimentation rate values in grid
  Nk=numk                              # number of precession frequency values in grid
  

#######################################################################################
#### (2) if check=T, ensure resolution requirements are satisfied
#######################################################################################
  if(check)
   {
    prior <- .getpriorpdfastropara(tGa,mvopt)                  # parameters for priors, excluding u
# check minimum and maximum sedimentation rates. sedmin is in meters/kyr, dx is in meters.
    NyqFreq=umin/(2*dx)
# freqHigh identifies the frequency of the shortest period in the target, given the 
#   stated uncertainties in g and k. we add ftol for the maximum possible upper half power 
#   point in the taner filtering.  we will take 2 standard deviations of mean as limit.
    eop <- .compeopfreq(c(prior$snmean+2*prior$sigma,0),mvopt)
    pvP=eop$pv*fconv
    freqHigh=max(pvP)+ftol
    if(freqHigh>NyqFreq)
     {
      umin = 2*dx*freqHigh
      if(verbose) cat("\n * WARNING: minimum sedimentation rate is too low.\n")
      if(verbose) cat("            -> sedmin reset to",100*umin,"cm/ka\n")
     }
# check maximum sedimentation rate. sedmax is in meters/kyr. dx is in meters.
    RayFreq = umax/(Ndv*dx)
# freqLow identifies the frequency of the longest period in the target, given the 
#   stated uncertainties in g. we will not worry about ftol here. we will take 2 standard 
#   deviations of mean as the nominal limit.
    eop <- .compeopfreq(c(prior$snmean-2*prior$sigma,0),mvopt)
    evP=eop$ev*fconv
    freqLow=min(evP)
    if(RayFreq>freqLow)
     {
      umax = Ndv*dx*freqLow
      if(verbose) cat("\n * WARNING: maximum sedimentation rate is too high.\n")
      if(verbose) cat("            -> sedmax reset to",100*umax,"cm/ka\n")
     }
    if(umin>umax) stop("**** ERROR: sedmin > sedmax. TERMINATING NOW!")
   }

#######################################################################################
#### (3) compute fft, set additional variables/constants, set-up pdfpara                
#######################################################################################
# calculate fourier coefficients for the stratigraphic series. the dataframe fc_dat
#  contains three columns: frequency, real fourier coeff, imaginary fourier coeff
  fc_dat <- periodogram(dat,padfac=2,demean=T,detrend=F,output=2,nrm=0,genplot=F,check=F,verbose=F)

# set additional variables/constants
# compute prior values of g's, s's, and k (arcsec/yr)
  if(!check) prior <- .getpriorpdfastropara(tGa,mvopt)         # parameters for priors, excluding u
# g's and s's are set to their prior mean value, which may differ from the
# mu parameter listed in Hoang et al. 2021
# here we are initializing mv to the prior means
  mv = prior$snmean                    # g's, s's, k (and u, to be added later)
  priormeanastro = prior$snmean 
  priorsigmaastro = prior$sigma
# add element for u (set to zero)
  mv=append(mv,0)

# compute bounds of precession frequency k if they are not given
  ik <- .getmvindex('k',mvopt)
  iu <- .getmvindex('u',mvopt)
  if (is.null(kmin) || is.null(kmax))
   {
    kdelta=min(6,4*priorsigmaastro[ik])                        # bound prior uncertainties in k. we set this to a minimum of 6, but could be removed or modified
    kmin=mv[ik]-kdelta
    kmax=mv[ik]+kdelta
   }

# grid of sedimentation rates (m/kyr) and precession frequencies
# (arcsec/yr) for evaluation of likelihood and posterior pdf
  uv=seq(umin,umax,length.out=Nu)
  kv=seq(kmin,kmax,length.out=Nk)
  
# prior pdfs of u and k for values in grid vectors uv and kv
  kpriorpdf <- .comppriorpdf1('k',kv,mvopt,tGa,umin,umax)
  upriorpdf <- .comppriorpdf1('u',100*uv,mvopt,tGa,100*umin,100*umax)  # for u in cm/kyr

# compute period of shortest eccentricity cycle and taudenv
  ev <- .compeopfreq(mv,mvopt)$ev
  periode=1/(fconv*max(ev))                                    # period (kyr) of shortest eccentricity cycle
  umid=(umin+umax)/2                                           # reference sedimentation rate (m/kyr)
  taudenv = umid*periode/(4*dx)                                # tau for envelope residuals at u=umid  (renamed from Matlab taudenvresid)
  cat('\n * taudenv:', taudenv,"\n")

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
  pdfpara = double(18)
  names(pdfpara) = c('tGa','P','roll','deltafp','nolikfdenv','mvopt','iu','umin','umax',
                      'ik','kmin','kmax','alpha','Nu','Nk','Nacf','taudenv','ar2')
  pdfpara['tGa'] = tGa                 # age in billions of years
  pdfpara['P'] = 2                     # order of AR(P) process fitted to the residuals
  pdfpara['roll'] = roll               # roll-off rate in Taner filter (dimensionless)
  pdfpara['deltafp'] = ftol            # default delta frequency for Taner filter boundaries

#  pdfpara['likfonly'] = 0             # if nonzero, use the likelihood as the target pdf (option not being used)
  if(env) pdfpara['nolikfdenv'] = 0    # if nonzero, omit the envelope likelihood
  if(!env) pdfpara['nolikfdenv'] = 1   # if nonzero, omit the envelope likelihood
  pdfpara['mvopt'] = mvopt             # variable that defines the astro. parameters in mv
  pdfpara['iu'] =  iu                  # index of u (sedimentation rate) in mv
  pdfpara['umin'] = umin               # minimum sedimentation rate to evaluate (m/kyr)
  pdfpara['umax'] = umax               # maximum sedimentation rate to evaluate (m/kyr)
  pdfpara['ik'] = ik                   # index of k (precession frequency) in mv
  pdfpara['kmin'] = kmin               # minimum k to evaluate (arcsec/yr)
  pdfpara['kmax'] = kmax               # maximum k to evaluate (arcsec/yr)
  pdfpara['alpha'] = 0.95              # probability for credible interval of u, k
  alpha = 0.95  # copy for present script

  pdfpara['Nu'] = numsed               # number of sedimentation rate values in grid
  pdfpara['Nk'] = numk                 # number of precession frequency values in grid
  pdfpara['Nacf'] = 21                 # number of lags in autocorrelation of AR(2) driving noise
  Nacf = 21     # copy for present script
  
  pdfpara['ar2'] = 4                   # option for ar.burg model fitting
  pdfpara['taudenv'] = taudenv         # tau for envelope residuals at u=umid
#   deltaz = dx                        # Matlab deltaz omitted, using dx from above

#######################################################################################
#### (3) analysis               
#######################################################################################

# ========================================================================
# ====== likelihoods for data and envelope, prior and posterior pdf ====== 
# ========================================================================

# set up progress display
  if(verbose)
   {
    cat("\n * Posterior pdf calculation progress:\n")
    cat("0%       25%       50%       75%       100%\n")
    progress <- utils::txtProgressBar(min = 0, max = Nk, style = 1, width=43)
   }
# initialize output matrices. 
# rows evaluate different precession frequencies, columns evaluate different sedimentation rates.
  loglikfdvimage = matrix(0,nrow=Nk,ncol=Nu)
  loglikfdenvimage = matrix(0,nrow=Nk,ncol=Nu)
  loglikfimage = matrix(0,nrow=Nk,ncol=Nu)
  logpriorpdfimage = matrix(0,nrow=Nk,ncol=Nu)
  logpostpdfimage = matrix(0,nrow=Nk,ncol=Nu)
  phi1image = matrix(0,nrow=Nk,ncol=Nu)
  phi2image = matrix(0,nrow=Nk,ncol=Nu)
#  maxloglikf = -realmax
  maxloglikf = -.Machine$double.xmax
  maxloglikfdv = -.Machine$double.xmax
  maxloglikfdenv = -.Machine$double.xmax
#  maxlogpostpdf = -realmax
  maxlogpostpdf = -.Machine$double.xmax

# loop over k, u
  for (i in 1:Nk)
   {
    if(verbose) utils::setTxtProgressBar(progress, i)   
    for (j in 1:Nu)
     {
      mv[ik]=kv[i]
      mv[iu]=uv[j]
# log-likelihood of data and envelope
#     dat contains dv,dz; pdfpara contains mvopt,iu,taudenv,P,roll,deltafp
      loglik <- .comploglikf(mv,dat,fc_dat,pdfpara)
      loglikfdvimage[i,j] = loglik$loglikfdv
      loglikfdenvimage[i,j] = loglik$loglikfdenv
      phiv=loglik$phivest
      phi1image[i,j]=phiv[1]
      phi2image[i,j]=phiv[2]
# log-total likelihood
      if (pdfpara['nolikfdenv']==1) loglikfimage[i,j]=loglikfdvimage[i,j] # ignore envelope likelihood
      if (pdfpara['nolikfdenv']==0) loglikfimage[i,j]=loglikfdvimage[i,j]+loglikfdenvimage[i,j]
# unnormalized log-prior pdf (k only, prior of u is uniform)
      logpriorpdfimage[i,j]=log(kpriorpdf[i])
# unnormalized log-posterior pdf
      logpostpdfimage[i,j]=logpriorpdfimage[i,j]+loglikfimage[i,j]
# track ML and MAP values of parameters
      if (loglikfdvimage[i,j]>maxloglikfdv)
       {
        jmaxlikdv=uv[j]
        imaxlikdv=kv[i]
        maxloglikfdv=loglikfdvimage[i,j]
       }
      if (loglikfdenvimage[i,j]>maxloglikfdenv)
       {
        jmaxlikdenv=uv[j]
        imaxlikdenv=kv[i]
        maxloglikfdenv=loglikfdenvimage[i,j]
       }
      if (loglikfimage[i,j]>maxloglikf)
       {
        uml=uv[j]
        kml=kv[i]
        phivml=phiv
        maxloglikf=loglikfimage[i,j]
       }
      if (logpostpdfimage[i,j]>maxlogpostpdf)
       {
        umap=uv[j]
        kmap=kv[i]
        phivmap=phiv
        maxlogpostpdf=logpostpdfimage[i,j]
       }
     }
   } 
# close progress display
  close(progress)     

# ===== compute rescaled likelihood and posterior pdf ===== 
# rescale log-likelihoods for data and data envelope by setting their
# mode to 0
  maxloglikf=max(loglikfdvimage)
  loglikfdvimage=loglikfdvimage-maxloglikf
  maxloglikf=max(loglikfdenvimage)
  loglikfdenvimage=loglikfdenvimage-maxloglikf

# rescale total log-likelihood (mode = 0) and total likelihood (mode = 1)
  maxloglikf=max(loglikfimage)
  loglikfimage=loglikfimage-maxloglikf
  likfimage=exp(loglikfimage)

# rescale log-posterior pdf (mode = 0) and posterior pdf (mode = 1)
  logpostpdfmap=max(logpostpdfimage)                           # MAP value of log-posterior pdf
  logpostpdfimage=logpostpdfimage-logpostpdfmap
  postpdfimage=exp(logpostpdfimage)

# ========================================================================
# ============================ output results ============================ 
# ========================================================================
# output prior pdf of k parameters
  if(verbose)
   {
    cat('\n========= Prior PDF =============================================\n')
    cat('  \tmean  ,\tst.dev.\n')
    cat(sprintf('k:\t%.5f  ,\t%.5f\n',priormeanastro[ik],priorsigmaastro[ik]))
   }
# output marg. lik. parameters for k
  kmargunnorm = rowSums(likfimage)                               # unnormalized
  kmarg <- .comppdfpara(kv,kmargunnorm,alpha)
  kmarglikmean = kmarg$xmean
  kmargliksdev = kmarg$xsdev
  kmargliklow = kmarg$xlow
  kmarglikhigh = kmarg$xhigh
# do the same for u in a column vector  
  umargunnorm = colSums(likfimage)
  umarg <- .comppdfpara(100*uv,umargunnorm,alpha)
  umarglikmean = umarg$xmean
  umargliksdev = umarg$xsdev
  umargliklow = umarg$xlow
  umarglikhigh = umarg$xhigh
  if(verbose)
   {
    cat('\n========= Marginal likelihood ===================================\n')
    cat(sprintf('  \tmax.lik.  ,\tmean  ,\tst.dev.  ,\t%g%% interval\n',100*alpha))
    cat(sprintf('u:\t%.5f  ,\t%.5f  ,\t%.5f  ,\t%.5f - %.5f\n',100*uml,umarglikmean,umargliksdev,umargliklow,umarglikhigh))
    cat(sprintf('k:\t%.5f  ,\t%.5f  ,\t%.5f  ,\t%.5f - %.5f\n',kml,kmarglikmean,kmargliksdev,kmargliklow,kmarglikhigh))
   }
# output marg. post. parameters
  kmargunnorm = rowSums(postpdfimage)                          # unnormalized
  kmarg <- .comppdfpara(kv,kmargunnorm,alpha)
  kmargpostmean = kmarg$xmean
  kmargpostsdev = kmarg$xsdev
  kmargpostlow = kmarg$xlow
  kmargposthigh = kmarg$xhigh
  kmargpostpdf = kmarg$pdfx
# do the same for u in a column vector
  umargunnorm=colSums(postpdfimage)
  umarg <- .comppdfpara(100*uv,umargunnorm,alpha)
  umargpostmean = umarg$xmean
  umargpostsdev = umarg$xsdev
  umargpostlow =  umarg$xlow
  umargposthigh = umarg$xhigh
  umargpostpdf = umarg$pdfx
  if(verbose)
   {
    cat('\n========= Marginal posterior PDF ================================\n')
    cat(sprintf('  \tMAP  ,\tmean  ,\tst.dev.  ,\t%g%% interval\n',100*alpha))
    cat(sprintf('u:\t%.5f  ,\t%.5f  ,\t%.5f  ,\t%.5f - %.5f\n',100*umap,umargpostmean,umargpostsdev,umargpostlow,umargposthigh))
    cat(sprintf('k:\t%.5f  ,\t%.5f  ,\t%.5f  ,\t%.5f - %.5f\n',kmap,kmargpostmean,kmargpostsdev,kmargpostlow,kmargposthigh))
   }
# ===== R^2 of data fit for MAP value of u, k =====
  mvmap=mv
  mvmap[ik]=kmap
  mvmap[iu]=umap
# Gdv for MAP value of mv
  tv=(dat[,1]-dat[1,1])/umap           # kyr
  Gdv <- .compGmats(mvmap,mvopt,tv)$Gdv
# compute R^2 values for data fit (no envelope) for MAP value of u, k
  sumsqrdv=sum(dat[,2]^2)
  GdvTGdv=t(Gdv)%*%Gdv
  afitv=solve(GdvTGdv,t(Gdv)%*% dat[,2])                       # least-squares fit coeffs.
  rsqdatamap <- .comprsq(dat[,2],sumsqrdv,Gdv,afitv)
# determine indices of columns of Gdv that correspond to eccentricity,
# obliquity, climatic precession and corresponding R^2 values
  eop <- .compeopfreq(mvmap,mvopt)
  ev= eop$ev
  ov= eop$ov
  pv= eop$pv
  nev=length(ev)
  je= 1:(nev*2)
  rsqdatamapecc <- .comprsq(dat[,2],sumsqrdv,Gdv,afitv,je)
  if(any(is.na(ov)))
   {
    nov=0
    rsqdatamapobl=0
   }else{ 
    nov=length(ov)
    jo=(nev*2)+(1:(nov*2))
    rsqdatamapobl <- .comprsq(dat[,2],sumsqrdv,Gdv,afitv,jo)
   }
  if (any(is.na(pv)))
   {
    npv=0
    rsqdatamappre=0
   }else{
    npv=length(pv)
    jp=((nev+nov)*2)+(1:(npv*2))
    rsqdatamappre <- .comprsq(dat[,2],sumsqrdv,Gdv,afitv,jp)
   }
  if (2*(nev+nov+npv) != dim(Gdv)[2]) stop('Mismatch in calculating indices of columns of Gdv') 
  if(verbose)
   {   
# report R^2 values
    cat('\n========= MAP R^2 values ========================================\n')
    cat(sprintf('all astronomical cycles:  %.3f\n',rsqdatamap))
    cat(sprintf('eccentricity only:        %.3f\n',rsqdatamapecc))
    if (nov>0) cat(sprintf('obliquity only:           %.3f\n',rsqdatamapobl))
    if (npv>0) cat(sprintf('climatic precession only: %.3f\n',rsqdatamappre))
    cat('                          -----\n')
    cat(sprintf('check sum (e+o+p):        %.3f\n', rsqdatamapecc+rsqdatamapobl+rsqdatamappre))
   }
# report values of phi1, phi2
  arpfit <- .arpestim(dat[,2],pdfpara)                         # phi1 and phi2 estimated on original data
  phivd = arpfit$phiv
  wv = arpfit$wv                                               # first P elements are NA
  wv = wv[(pdfpara['P']+1):Ndv]

  phi1mean=mean(phi1image)                                     # mean phi1, phi2 of residuals
  phi2mean=mean(phi2image)
  if(verbose)
   {
    cat('\n========= AR(2) process coefficients ============================\n')
    cat('                       \tphi1  ,\tphi2\n')
    cat(sprintf('observed data:         \t%.2f  ,\t%.2f\n',phivd[1],phivd[2]))
    cat(sprintf('mean of all residuals: \t%.2f  ,\t%.2f\n',phi1mean,phi2mean))
    cat(sprintf('residuals at ML point: \t%.2f  ,\t%.2f\n',phivml[1],phivml[2]))
    cat(sprintf('residuals at MAP point:\t%.2f  ,\t%.2f\n',phivmap[1],phivmap[2]))
   }
# use phi1, phi2 at MAP point
  phi1=phivmap[1]
  phi2=phivmap[2]

# autocorrelation of driving noise of AR(2) process fitted to the data dv
#  r=xcorr(wv,Nacf,'normalized');
  acfwv=acf(wv,lag.max=Nacf,type="correlation",demean=TRUE,plot=FALSE)$acf
#acfwvbound=1.96/sqrt(Ndv); % 95% bound on acf of white noise
  acfwvbound=1.96/sqrt(Ndv-2)                                  # 95% bound on acf of white noise

# data periodogram and AR(2) spectrum
  datT=data.frame(cbind(tv,dat[,2]))
#[fv,pgdv]=pgram(dv,deltat,zeromult); % periodogram of observed data
  pgdat=periodogram(datT,padfac=2,demean=T,detrend=F,output=1,nrm=1,genplot=F,check=F,verbose=F)
  fv=pgdat[,1]
  pgdv=pgdat[,3]
  Nfv=length(fv)
# spectrum of AR(2) process (p. 123 of Chatfield 1975)
  omega=pi*fv/fv[Nfv]                                          # angular frequencies for deltat=1
  cosomega=cos(omega)
  cos2omega=cos(2*omega)
  ar2sp=(1+phi1^2+phi2^2-2*phi1*(1-phi2)*cosomega-2*phi2*cos2omega)^(-1)
  norm=sum(pgdv)/sum(ar2sp)
  ar2sp=norm*ar2sp

# Monte Carlo simulation (Nsim>0) is in separate function 'timeOptBSim'


# ========================================================================
# ============================      plot      ============================
# ========================================================================
  if(genplot)
   {
# set color palette
    ncolors=100
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

# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Figure 1: Log-likelihood and likelihood images for spectral and 
#           envelope fit (4 plots)
    dev.new(height = 5.9, width = 7.3, units = "in")
    layoutA <- layout(matrix(c(1,2,3,4), 2, 2, nrow=2)) 
#          bottom, left, top, right
    par(mar = c(3.2, 3.6, 2.1, 1))
# spectral log-likelihood image
    image.plot(uv*100,kv,t(loglikfdvimage),xlab="",ylab="",main="Log-likelihood of data",col=colPalette)
    mtext(c("Sedimentation rate (cm/kyr)"),side=1,line=2)
    mtext("Precession frequency (arcsec/yr)",side=2,line=2)
#    points(jmaxlikdv*100,imaxlikdv,pch=3,col="red",cex=1.8)
    points(jmaxlikdv*100,imaxlikdv,pch=19,col="red",cex=1.1)
    points(jmaxlikdv*100,imaxlikdv,pch=21,col="white",cex=1.3)
# envelope log-likelihood image
    image.plot(uv*100,kv,t(loglikfdenvimage),xlab="",ylab="",main="Log-likelihood of data envelope",col=colPalette)
    mtext(c("Sedimentation rate (cm/kyr)"),side=1,line=2)
    mtext("Precession frequency (arcsec/yr)",side=2,line=2)
#    points(jmaxlikdenv*100,imaxlikdenv,pch=3,col="red",cex=1.8)
    points(jmaxlikdenv*100,imaxlikdenv,pch=19,col="red",cex=1.1)
    points(jmaxlikdenv*100,imaxlikdenv,pch=21,col="white",cex=1.3)
# cross out envelope result if not used
    if (pdfpara['nolikfdenv'] == 1) 
     {
      lines(c(100*uv[1],100*uv[Nu]),c(kv[1],kv[Nk]),col="#BEBEBE5A",lwd=10)
      lines(c(100*uv[1],100*uv[Nu]),c(kv[Nk],kv[1]),col="#BEBEBE5A",lwd=10)
     }
# total log-likelihood image
    image.plot(uv*100,kv,t(loglikfimage),xlab="",ylab="",main="Total log-likelihood",col=colPalette)
    mtext(c("Sedimentation rate (cm/kyr)"),side=1,line=2)
    mtext("Precession frequency (arcsec/yr)",side=2,line=2)
#    points(uml*100,kml,pch=3,col="red",cex=1.8)    
    points(uml*100,kml,pch=19,col="red",cex=1.1) 
    points(uml*100,kml,pch=21,col="white",cex=1.3)
# total likelihood image
    image.plot(uv*100,kv,t(likfimage),xlab="",ylab="",main="Total likelihood",col=colPalette)
    mtext(c("Sedimentation rate (cm/kyr)"),side=1,line=2)
    mtext("Precession frequency (arcsec/yr)",side=2,line=2)
#    points(uml*100,kml,pch=3,col="red",cex=1.8) 
    points(uml*100,kml,pch=19,col="red",cex=1.1)
    points(uml*100,kml,pch=21,col="white",cex=1.3)  
    abline(h=k0,col="orange",lwd=2,lty=5)
    text(x=100*umid,y=k0,c("Present-day k"),col="white",font=2)
  
# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Figure 2: Log-posterior and posterior PDF images, marginal posterior 
#           PDFs of u and k
    dev.new(height = 5.9, width = 7.3, units = "in")
    layoutA <- layout(matrix(c(1,2,3,4), 2, 2, nrow=2)) 
#          bottom, left, top, right
    par(mar = c(3.2, 3.6, 2.1, 1))
 
# log-posterior pdf image
    image.plot(uv*100,kv,t(logpostpdfimage),xlab="",ylab="",main="Log-posterior PDF",col=colPalette)
#    points(umap*100,kmap,pch=3,col="red",cex=1.8)    
    points(umap*100,kmap,pch=19,col="red",cex=1.1)    
    points(umap*100,kmap,pch=21,col="white",cex=1.3)  
    mtext(c("Sedimentation rate (cm/kyr)"),side=1,line=2)
    mtext("Precession frequency (arcsec/yr)",side=2,line=2)
# posterior pdf image
    image.plot(uv*100,kv,t(postpdfimage),xlab="",ylab="",main="Posterior PDF",col=colPalette)
    points(umap*100,kmap,pch=19,col="red",cex=1.1) 
    points(umap*100,kmap,pch=21,col="white",cex=1.3) 
    mtext(c("Sedimentation rate (cm/kyr)"),side=1,line=2)
    mtext("Precession frequency (arcsec/yr)",side=2,line=2)
    abline(h=k0,col="orange",lwd=2,lty=5)
    text(x=100*umid,y=k0,c("Present-day k"),col="white",font=2)
    
# x-coordinate limits for marginal posterior plots (make sure plot min. is
# >= xmin and plot max. is <= xmax)
    uvplotmin=max(100*umin,umargpostmean-4*umargpostsdev)
    uvplotmax=min(100*umax,umargpostmean+7*umargpostsdev)      # space for legend
    kvplotmin=max(kmin,kmargpostmean-4*kmargpostsdev)
    kvplotmax=min(kmax,kmargpostmean+4*kmargpostsdev)
# y-coordinate limits for marginal posterior plots
    uvplotmax2=max(umargpostpdf,upriorpdf)
    kvplotmax2=max(kmargpostpdf,kpriorpdf)

# marginal posterior pdf of u
    plot(100*uv,umargpostpdf,type="n",xlab="",ylab="",main="Marginal posterior PDF of sed. rate",ylim=c(0,uvplotmax2),col="red")
    polygon( c(100*uv,rev(100*uv)) , c(umargpostpdf,rep(0,Nu)) , col="#FF000050",border=NA )
    lines(100*uv,upriorpdf,col="blue",lwd=2)
    abline(v=umargpostmean,col="red",lty=2,lwd=2)
    points(100*umap,max(umargpostpdf),pch=19,col="red",cex=1.1)
    points(100*umap,max(umargpostpdf),pch=21,col="white",cex=1.3) 
    mtext(c("Sedimentation rate (cm/kyr)"),side=1,line=2)
    if(abs(100*umax-umargpostmean) <  abs(100*umin-umargpostmean)) 
     {
      legend(x="topleft",legend=c('Posterior PDF','MAP value','Posterior mean','Prior PDF'),col=c("#FF000050","red","red","blue"),lty=c(1,NA,2,1),lwd=c(8,NA,2,2),pch=c(NA,19,NA,NA),cex=0.7)
     }
    if(abs(100*umax-umargpostmean) >  abs(100*umin-umargpostmean)) 
     {
      legend(x="topright",legend=c('Posterior PDF','MAP value','Posterior mean','Prior PDF'),col=c("#FF000050","red","red","blue"),lty=c(1,NA,2,1),lwd=c(8,NA,2,2),pch=c(NA,19,NA,NA),cex=0.7)
     }
# check for multimodality of marginal posterior pdf of u
    numpeak <- length(peak(umargpostpdf,level=0.05*(max(umargpostpdf)),genplot=F,verbose=F)[,1]) 

# marginal posterior pdf of k
    plot(kv,kmargpostpdf,type="n",xlab="",ylab="",main="Marginal posterior PDF of k",ylim=c(0,kvplotmax2),col="red")
    polygon( c(kv,rev(kv)) , c(kmargpostpdf,rep(0,Nk)) , col="#FF000050",border=NA)
    lines(kv,kpriorpdf,col="blue",lwd=2)
    abline(v=kmargpostmean,col="red",lty=2,lwd=2)
    points(kmap,max(kmargpostpdf),pch=19,col="red",cex=1.1)
    points(kmap,max(kmargpostpdf),pch=21,col="white",cex=1.3)
    mtext("Precession frequency (arcsec/yr)",side=1,line=2)    

# ensure that tails of posteriors approach zero
   if (umargpostpdf[1] > 0.05*max(umargpostpdf) || umargpostpdf[Nu] > 0.05*max(umargpostpdf))   
    {
      cat("\n * WARNING: marginal posterior pdf for sedimentation rate may not be fully captured.\n")  
      cat("            Rerun analysis with a broader range for sedimentation rate.\n\n") 
    }
   if (kmargpostpdf[1] > 0.05*max(kmargpostpdf) || kmargpostpdf[Nk] > 0.05*max(kmargpostpdf))
    {
      cat("\n * WARNING: marginal posterior pdf for k may not be fully captured.\n")  
      cat("            Rerun analysis with a broader range for k.\n\n") 
    }

# if marginal posterior pdf of u is unimodal, conduct another check for multimodality of k
    if(numpeak == 1) numpeak <- length(peak(kmargpostpdf,level=0.05*(max(kmargpostpdf)),genplot=F,verbose=F)[,1]) 
    if (numpeak > 1)
     {
      cat("\n * WARNING: possible multimodal marginal posterior pdf detected.\n")  
      cat("            This may be an artifact of the analysis grid.\n")
      cat("            Rerun analysis with a finer grid for sedimentation rate and k.\n\n") 
     }      

# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# Figure 3: fit to the data for MAP value of u and k
   .plotmvdatafit(dat,fc_dat,mvmap,mvopt,iu,roll,pdfpara['deltafp'],pdfpara['nolikfdenv'])

# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Figure 4: AR(2) process fitted to the data
    dev.new(height = 7.3, width = 6.5, units = "in")
    layoutA <- layout(matrix(c(1,2), 2, 1)) 
#          bottom, left, top, right
    par(mar = c(3.2, 3.6, 2.1, 1))
  
    mv[ik]=kmax
    pv <- .compeopfreq(mv,mvopt)$pv
    astrofmax=fconv*max(pv)            # max. astronomical frequency considered
    fmax=min(3*astrofmax,tail(fv,1))                           # max. frequency to plot
    ifmax=which(fv>=fmax)[1]
    ifrange=1:ifmax

# note that there is no zero freq returned from periodogram
    plot(fv[ifrange],pgdv[ifrange],type="l",log="y",xlab="",ylab="",main="AR(2) process fit to the data")
    lines(fv[ifrange],ar2sp[ifrange],col="purple",lwd=2,lty=5)
    abline(v=astrofmax,lty=3,lwd=1.5)
    mtext(c("Frequency (cycles/kyr)"),side=1,line=2)
    mtext(c("Power"),side=2,line=2)
    legend(x="topright",legend=c('Power spectrum of data','Spectrum of AR(2) process','Max. astronomical frequency'),col=c("black","purple","black"),lty=c(1,5,3),lwd=c(1,2,1,5),cex=0.8,bty="n")

# ACF of AR(2) process driving noise                                                                   ! 
    lagv=0:Nacf                        # MATLAB extracts lags from 0-20, while R 0-21
    plot(lagv,acfwv,type="h",xlab="",ylab="",main="")
    points(lagv,acfwv,pch=19)
    abline(h=0,lty=5)
    mtext(c("Lag"),side=1,line=2)
    mtext("Autocorrelation",side=2,line=2)
    rect(-1,-acfwvbound,Nacf+1,acfwvbound,col="#BEBEBE50",border=NA)
    legend(x="topright",legend=c('95% bounds for white noise','Driving noise of AR(2) process'),col=c("#BEBEBE50","black"),lty=c(1,5,3),lwd=c(8,1),pch=c(NA,19),cex=0.8,bty="n")

# end genplot
  }
  
# output parameters for timeOptBSim
  out = double(17)
  names(out) = c('Ndv','dx','mvopt','tGa','umin','umax','Nu','kmin','kmax','Nk','phi1','phi2',
                    'rsqdatamap','rsqdatamapecc','rsqdatamapobl','rsqdatamappre','detrend')
  out['Ndv'] = Ndv
  out['dx'] = dx
  out['mvopt'] = mvopt
  out['tGa'] = tGa
  out['umin'] = umin
  out['umax'] = umax
  out['Nu'] = Nu
  out['kmin'] = kmin
  out['kmax'] = kmax
  out['Nk'] = Nk
  out['phi1'] = phi1
  out['phi2'] = phi2
  out['rsqdatamap'] = rsqdatamap
  out['rsqdatamapecc'] = rsqdatamapecc
  out['rsqdatamapobl'] = rsqdatamapobl
  out['rsqdatamappre'] = rsqdatamappre
  out['detrend'] = detrend
  return(invisible(out))  
 
# end function timeOptB
 }
 
