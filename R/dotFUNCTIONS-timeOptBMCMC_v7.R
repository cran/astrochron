### These functions are components of astrochron: An R Package for Astrochronology
### Copyright (C) 2026 Stephen R. Meyers
###
#######################################################################################
#### Definition of key .FUNCTIONS for timeOptBMCMC.R (SRM: Feb. 21, 2026)
#### This includes:
####   .getmvlabels, .getmvindex, .compeopfreq, .mcmcadapt, .mcmclogpdf
####   .getpriorpdfastropara, .complogpriorpdf, .skewnormalpdf, .arpestim
####   .comploglikf, .compdpred, .compGmats, .comppriorpdf1, .comppdfpara
####   .comppvalue, .comprsq, .plotmvdatafit
#######################################################################################
# ---------- FUNCTION getmvlabels ----------
# getmvlabels is an R translation of Alberto Malinverno's MATLAB function.
#  return all label(s) of al model vector parameters, depending on mvopt
#
# Return labels depending on mvopt. Assumptions:
# - The first 5 elements of mv are always g_1 to g_5
# - The last 2 elements of mv are always k and u, respectively
.getmvlabels <- function(mvopt)
 {
# full list of model parameter labels
#  mvlabelstr=c('g_1', 'g_2', 'g_3', 'g_4', 'g_5', 's_1', 's_2', 's_3',
#    's_4', 's_6', 'k', 'u')
  mvlabelstr=c('g1', 'g2', 'g3', 'g4', 'g5', 's1', 's2', 's3',
    's4', 's6', 'k', 'u')
# depending on mvopt, delete labels from full list
  if(mvopt==0) idel=0
  if(mvopt==1) idel=c(6, 7, 10)        # delete s_i entries except for s_3 and s_4
  if(mvopt==2) idel=6:10               # delete all s_i
  if (mvopt != 0) mvlabelstr = mvlabelstr[-idel]
  return(mvlabelstr)
 }   

# ---------- FUNCTION getmvindex ----------
# getmvindex is an R translation of Alberto Malinverno's MATLAB function.
# return index(es) corresponding to label(s) of model vector parameters.
# The value of mvlabel can be
# - 'g' or 's': return the indexes of the corresponding g_i and s_i
#   frequencies in the model vector
# - A single parameter, e.g., 'g_3', 's_4', 'k', 'u'
# - Multiple parameters, e.g., mvlabel={'g_1', 'g_2', 'g_3, ...}
#
# Return labels depending on mvopt. Assumptions:
# - The first 5 elements of mv are always g_1 to g_5
# - The last 2 elements of mv are always k and u, respectively
.getmvindex <- function(mvlabel,mvopt)
 {
# get list of model parameter labels given mvopt
  mvlabelstr <- .getmvlabels(mvopt)
# return index(es) of name(s) in mvlabelstr
  if (is.character(mvlabel) && length(mvlabel)==1 ) # if mvlabel has one entry
   {
    if(mvlabel == 'g') mvindex=1:5     # all g_i freqs, assumed to be always be elements 1-5
    if(mvlabel == 's')                
     {
      if(mvopt == 0) mvindex=6:10      # include all s_i freqs
      if(mvopt == 1) mvindex=6:7       # only include s_3 and s_4
      if(mvopt == 2) mvindex=NA        # no 's' freqs. in model vector
     }
    if(mvlabel != 'g' && mvlabel != 's') mvindex=which(mvlabel==mvlabelstr) # single parameter, e.g., 'g_3', 's_4', 'k', 'u'

   }else{                              # if mvlabel contains multiple entries
    nlabels=length(mvlabel)
    mvindex=double(nlabels)
    for (i in 1:nlabels)
     {
      thismvindex=which(mvlabel[i]==mvlabelstr)
      mvindex[i]=thismvindex
     }
   }
  return(mvindex)
 } 

# ---------- FUNCTION compeopfreq ----------
# compeopfreq is an R translation of Alberto Malinverno's MATLAB function.
# compute frequencies of eccentricity (in vector ev), obliquity (vector
# ov), and climatic precession (vector pv) for given fundamental Solar
# system frequencies (vectors gv and sv) and precession frequency (k),
# which are in the first 11 elements of the model parameter vector mv.
# Returned frequencies ev, ov, and pv have the same the units as gv, sv,
# and k (e.g., arcsec/yr or cycles/kyr).
#
# Modified to allow for different choices of secular
# frequencies g_i and s_i, indicated by mvopt. If there are no frequencies
# of a given type, the corresponding vector is set to an empty vector
# (length 0). Assumptions:
# - The first 5 elements of mv are always g_1 to g_5
# - The last 2 elements of mv are always k and u, respectively
# - The first 10 columns of Gdv contain sine and cosine terms for 5
#   eccentricity freqs. that are differences g_i - g_j
.compeopfreq <- function(mv,mvopt,label=FALSE)
 {
# determine gv (assumed to be always in the first 5 elements of mv) and
# corresponding eccentricity freqs. ev
  if (length(mv)<5) stop('Vector mv must have at least 5 g_i frequencies')
  gv=mv[1:5]
# 5 eccentricity frequencies g(i)-g(j), ordered from low to high
  ev=double(5)
  ev[1]=gv[2]-gv[5]
  ev[2]=gv[3]-gv[2]
  ev[3]=gv[4]-gv[2]
  ev[4]=gv[3]-gv[5]
  ev[5]=gv[4]-gv[5]
# obliquity frequencies
# depending on mvopt, extract sv,k from mv and determine
# ov (which may be empty and set to NA)
  if (mvopt == 0)                      # 5 obliquity and 5 climatic precession freqs.
   {
    if (length(mv) != 12) stop('Vector mv must have length 12')
    sv=mv[6:10]
    k=mv[11]        
    ov=double(5)                       # 5 obliquity frequencies ov=sv+k, ordered from low to high
    ov[1]=sv[5]+k                      # sv(5) is s_6 (s5=0 by definition)
    ov[2]=sv[3]+k
    ov[3]=sv[4]+k
    ov[4]=sv[2]+k
    ov[5]=sv[1]+k
   }
  if (mvopt == 1)                      # 2 obliquity and 5 climatic precession freqs.
   {
    if (length(mv) != 9) stop('Vector mv must have length 9')
    sv=mv[6:7]
    k=mv[8]
    ov=double(2)                       # 2 obliquity frequencies ov=sv+k
    ov[1]=sv[1]+k                      # s_3 + k
    ov[2]=sv[2]+k                      # s_4 + k
   }
  if (mvopt == 2)                      # No obliquity and 5 climatic precession freqs.
   {
    if (length(mv) != 7) stop('Vector mv must have length 9')
    k=mv[6]
    ov=NA                              # No obliquity frequencies, ov set to an empty vector
   }
# five climatic precession frequencies are always g_i + k
  pv=double(5)
  pv[1]=gv[5]+k                        # ordered from low to high
  pv[2]=gv[1]+k
  pv[3]=gv[2]+k
  pv[4]=gv[3]+k
  pv[5]=gv[4]+k

  if(!label) return(list(ev=ev,ov=ov,pv=pv))
# if requested as output, return frequency labels
  if(label)
   {
     evlabel=c('g2-g5','g3-g2','g4-g2','g3-g5','g4-g5')
     pvlabel=c('g5+k','g1+k','g2+k','g3+k','g4+k')
     if (mvopt == 0) ovlabel=c('s6+k','s3+k','s4+k','s2+k','s1+k')
     if (mvopt == 1) ovlabel=c('s3+k','s4+k')
     if (mvopt == 2) ovlabel=NA
     return(list(ev=ev,ov=ov,pv=pv,evlabel=evlabel,ovlabel=ovlabel,pvlabel=pvlabel))
   }
 }

# ---------- FUNCTION mcmcadapt ----------
# mcmcadapt is an R translation of Alberto Malinverno's MATLAB function.
# It is the generic Metropolis-within-Gibbs adaptive algorithm, as described in
#   Section 3 of:
#     Roberts, G.O., Rosenthal, J.S., 2009. Examples of Adaptive MCMC. Journal
#     of Computational and Graphical Statistics 18, 349-367.
#     https://doi.org/10.1198/jcgs.2009.06134
# mcmcadapt calls the function mcmclogpdf(mv,dv,pdfpara), which returns the value of the
#   target PDF, where pdfpara is a structure that contains parameters necessary for
#   the calculation of the target pdf. These parameters are specific to each
#   problem and must be provided as an argument to mcmcadapt()
# Required input arguments:
# - mv: Starting value of model parameter vector
# - dat (data frame): 
#   - depth 'vector' (pdfpara.zv of MATLAB)
#   - data 'vector'  (dv of MATLAB), detrended and normalized stratigraphic data
# - fc_dat (data frame)
#   - 'vector' of spatial frequencies
#   - 'vector' of real fourier coeff of dat (pdfpara.ftdv of MATLAB)
#   - 'vector' of imaginary fourier coeff of dat (pdfpara.ftdv of MATLAB)
# - prior: list of parameters defining the priors (excluding u), derived from getpriorpdfastropara 
#    (mu,sigma,alpha,a1,mu1,sigma1,snmean)
# - pdfpara: Structure with parameters for calculation of target pdf
# Optional input arguments (see default values in code):
# - sigmamvcand: Starting value of standard deviations of candidate model
#   vector perturbations
# - niter: Total number of iterations
# - niterinbatch: Number of iterations in each batch. At the end of each
#   batch, sigmamvcand is updated
# - savefile: save output file (T or F)
# Output is in 3 matrices:
# - mcmcadapt_mv: Sampled values of parameters
# - mcmcadapt_pacc: Probability of acceptance
# - mcmcadapt_sigmamvcand: Standard deviations of candidate model
#   vector perturbations
.mcmcadapt <- function (mv,dat,fc_dat=NULL,prior=NULL,pdfpara=NULL,sigmamvcand=NULL,niter=10000,niterinbatch=50,savefile=FALSE)
 {
  M = length(mv)		               # number of model parameters
# check input and compute number of batches
  if (niter %% niterinbatch !=0)  stop('mcmcadapt: niter must be an integer multiple of niterinbatch')
  nbatches = niter/niterinbatch
# for the case when sigmamvcand=NULL, set sigmamvcand values to unity
  if (is.null(sigmamvcand)) sigmamvcand = matrix(1,nrow=M,ncol=1)
# initialize output matrices, one row per iteration
  mcmcadapt_mv = matrix(NA,nrow=niter,ncol=M+1)                # add column for logpdf
  colnames(mcmcadapt_mv) <- c('logpdf',.getmvlabels(pdfpara['mvopt']))
  mcmcadapt_pacc = matrix(NA,nrow=nbatches,ncol=M)
  colnames(mcmcadapt_pacc) <- .getmvlabels(pdfpara['mvopt'])
  mcmcadapt_sigmamvcand = matrix(NA,nrow=nbatches,ncol=M)
  colnames(mcmcadapt_sigmamvcand) <- .getmvlabels(pdfpara['mvopt'])
# compute the initial value of log-target pdf
  logpdf <- .mcmclogpdf(mv,dat,fc_dat,prior,pdfpara)
  iter = 0		                       # iteration counter
# output starting mv and logpdf to output file with sampled model parameter vectors
  if(savefile) write.table(file="mcmc-mv.csv", matrix(data=c(iter, logpdf, mv),nrow=1,ncol=14), sep = ",", row.names = FALSE, col.names= FALSE, append=TRUE) 
# set up progress display
  cat("\n mcmcadapt progress:\n")
  cat("\n0%       25%       50%       75%       100%\n")
  progress <- utils::txtProgressBar(min = 0, max = nbatches, style = 1, width=43)
# begin main loop
  for (ibatch in 1:nbatches)
   {
    utils::setTxtProgressBar(progress, ibatch)
    nacc = matrix(0,nrow=M,ncol=1)	   # initialize number of accepted candidates in the batch
# begin batch loop
    for (iterinbatch in 1:niterinbatch)
     {
# parameter loop - try updating each parameter
      mvindex <- sample(1:M)	       # randomized indices from 1 to M
      for (j in 1:M)
       {
        i = mvindex[j]	               # select i-th model parameter in random order
        mvcand = mv
        mvcand[i] = mvcand[i]+sigmamvcand[i]*rnorm(1)	       # perturb i-th parameter
        logpdfcand <- .mcmclogpdf(mvcand,dat,fc_dat,prior,pdfpara)
        alpha = exp(logpdfcand-logpdf)					       # Metropolis acceptance prob.
        if (alpha>1 || runif(1)<alpha)					       # candidate accepted
         {
          mv = mvcand
          logpdf = logpdfcand
          nacc[i] = nacc[i]+1
         }
       }
# increment iteration counter
      iter = iter+1
# save values of logpdf, mv for this iteration (note that this excludes iter=0, unlike mcmc-mv.csv)
      mcmcadapt_mv[iter,1] = logpdf
      mcmcadapt_mv[iter,2:(M+1)] = mv[1:M]
# output iter, logpdf, mv to file mcmc-mv.csv
      if(savefile) write.table(file="mcmc-mv.csv", matrix(data=c(iter, logpdf, mv),nrow=1,ncol=14), sep = ",", row.names = FALSE, col.names= FALSE, append=TRUE)   
# end batch loop
     }      
    pacc = nacc/niterinbatch	       # probability of i-th parameter being accepted
# save values of pacc for this iteration, which is the end of the batch
    mcmcadapt_pacc[ibatch,1:M] = pacc[1:M]
# output iter, pacc to file mcmc-pacc.csv
    if(savefile) write.table(file="mcmc-pacc.csv", matrix(data=c(iter, pacc),nrow=1,ncol=13), sep = ",", row.names = FALSE, col.names= FALSE, append=TRUE)
# update sigmamvcand
    ls = log(sigmamvcand)
    delta = min(0.01,1/sqrt(iter))
    for (i in 1:M)
     {
      if (pacc[i]<0.44) 
       {
        ls[i] = ls[i]-delta	           # decrease ls(i) by delta
       } else {
        ls[i] = ls[i]+delta	           # increase ls(i)
       }
     }
    sigmamvcand = exp(ls)  
# save new updated values of sigmamvcand
    mcmcadapt_sigmamvcand[ibatch,1:M] = sigmamvcand[1:M]
# output iter, sigmamvcand to file mcmc-sigmamvcand.csv    
    if(savefile) write.table(file="mcmc-sigmamvcand.csv", matrix(data=c(iter, sigmamvcand),nrow=1,ncol=13), sep = ",", row.names = FALSE, col.names= FALSE, append=TRUE)   
# end main loop
   }
# close progress display
  close(progress)
# return results
  return(list(mv=mcmcadapt_mv,pacc=mcmcadapt_pacc,sigmamvcand=mcmcadapt_sigmamvcand))
# end function mcmcadapt
 }

# ---------- FUNCTION mcmclogpdf ----------
# mcmclogpdf is an R translation of Alberto Malinverno's MATLAB function.
# it returns the log-target pdf that is used in mcmcadapt.
# Input arguments:
# - mv: Model parameter vector
# - dat (data frame): 
#   - depth 'vector' (pdfpara.zv of MATLAB)
#   - data 'vector'  (dv of MATLAB), detrended and normalized stratigraphic data
# - fc_dat (data frame)
#   - 'vector' of spatial frequencies
#   - 'vector' of real fourier coeff of dat (pdfpara.ftdv of MATLAB)
#   - 'vector' of imaginary fourier coeff of dat (pdfpara.ftdv of MATLAB)
# - prior: list of parameters defining the priors (excluding u), derived from getpriorpdfastropara 
#    (mu,sigma,alpha,a1,mu1,sigma1,snmean)
# - pdfpara: Structure with parameters for calculation of target pdf.
#   The contents of this structure vary depending on the kind of pdf 
#   that is being sampled to allow for different applications.
.mcmclogpdf <- function(mv,dat,fc_dat,prior,pdfpara)
 {
  logpriorpdf <- .complogpriorpdf(mv,prior,pdfpara['mvopt'],pdfpara['tGa'],pdfpara['umin'],pdfpara['umax'],pdfpara['iu'],pdfpara['ik'])
  if (pdfpara['prioronly'] == 1)
   {
    logpdf = logpriorpdf
   } else {
    loglik <- .comploglikf(mv,dat,fc_dat,pdfpara)
    loglikfdv = loglik$loglikfdv
    loglikfdenv = loglik$loglikfdenv
    logpdf = loglikfdv+loglikfdenv+logpriorpdf                 # default
    if (pdfpara['nolikfdenv'] == 1)
     {
      if (pdfpara['likfonly'] == 1)
       {
        logpdf <- loglikfdv
       } else {
        logpdf <- loglikfdv+logpriorpdf
       }
     }
    if (pdfpara['likfonly'] == 1)
     {
      logpdf <- loglikfdv+loglikfdenv
     }
   } 
  return(logpdf)
 }

# ---------- FUNCTION getpriorpdfastropara ----------
# getpriorpdfastropara is an R translation of Alberto Malinverno's MATLAB function.
# compute parameters of prior pdf of astronomical parameters (5 g_i's,
# 5 s_i's, and k; 11 parameters in total) at age tGa. The parameters are
# those of a skew normal pdf (mu, sigma, and alpha) plus the parameters of
# a normal pdf (a1, mu1, sigma1) that define a secondary peak in the pdfs
# of g_4 and s_3. If alpha=0, the skew normal is a normal pdf. 
# NOTE: This function does not return parameters of the prior pdf 
# of u, which is a uniform distribution.
#
# Modified to only return prior parameters for some secular
# frequencies g_i and s_i, as indicated by mvopt. Assumptions:
# - The first 5 elements of mv are always g_1 to g_5
# - The last 2 elements of mv are always k and u, respectively
.getpriorpdfastropara <- function(tGa,mvopt)
 {
# ===== parameters of prior pdf of g_i and s_i (5 + 5) =====
# values in Table 2 of Hoang et al. 2021
# NOTE: mu_0 for s_6 is erroneously set to -2.634787 in Table 2
  amu=c(5.759, 7.448, 17.269, 17.896, 4.257454,                # g_1 to g_5 mu_0 in Table 2
    -5.652, -6.709, -18.773, -17.707, -26.34787)               # s_1 to s_6 mu_0 (no s_5)
  bmu=c(0.006, -0.004, 0.002, 0.005, -2.1E-6,                  # time-dependent change of mu_0
    -0.032, 0.030, 0.009, 0.013, 1.5E-5)
  asigmasq=c(3.37E-2, 4.17E-4, 6.63E-3, 6.88E-3, 4.63E-10,     # a in Table 2
    2.68E-2, 1.20E-1, 2.86E-2, 1.19E-2, 1.21E-8)
  bsigmasq=c(0.52, 0.70, 0.43, 0.41, 0.88,                     # exponent b of time-dependent sigmasq
    0.83, 0.76, 0.56, 0.68, 0.85)
  aalpha=c(-2.25, 1.38, 0, 0, 0,                               # skewness parameter alpha_0 in Table 2
    1.12, -2.94, -3.40, -1.73, 0)
  balpha=c(-0.5, 0.21, 0, 0, 0,                                # time-dependent change of alpha_0
   0.16, -1.23, -0.08, -0.28, 0)
# ===== parameters of prior pdf of k =====
# coeffs. of polynomial fit to prior mean of k, 0-3.3 Ga
# (from Figure 6 of Farhat et al. 2022, listed in AstroGeo)
  pmeank=c(2.4322, -11.2346, 13.0658, 23.1305, 50.4677)
# coeffs. of polynomial fit to multiplier of mean that gives the
# standard dev. of k, 0-3.3 Ga (from Waltham 2015 calculator)
  psigmamultk=c(0.0030, -0.0262, 0.0962, 0)
# ===== allocate output vectors =====
  nfreq=11
  mu=double(nfreq)
  sigma=double(nfreq)
  alpha=double(nfreq)
  a1=double(nfreq)
  mu1=double(nfreq)
  sigma1=double(nfreq)
  snmean=double(nfreq)
# ===== Solar system frequencies g_i and s_i (5 + 5) ====
  for(i in 1:10)
   {
    mu[i]=amu[i]+bmu[i]*tGa
    sigmasq=asigmasq[i]*tGa^(bsigmasq[i])
    sigma[i]=sqrt(sigmasq)
    alpha[i]=aalpha[i]+balpha[i]*tGa                           # skewness parameter of skew normal pdf
    delta=alpha[i]/sqrt(1+alpha[i]^2)
    snmean[i]=mu[i]+sigma[i]*delta*sqrt(2/pi)                  # mean of skew-normal pdf
    if (i==4)                                                  # secondary mode of g_4
     {
      a1[i]=0.11-0.012*tGa
      mu1[i]=17.6755
      sigma1[i]=sqrt(0.0034)
     } 
    if (i==8)                                                  # secondary mode of s_3
     {
      a1[i]=0.023
      mu1[i]=-18.5256
      sigma1[i]=sqrt(0.0028)
     }   
   }
# ===== precession frequency k ======
  i=nfreq
# mu(i)=polyval(pmeank,tGa);
# this inspired by polyval in pracma
  mu[i]= outer(tGa, 4:0, "^") %*% pmeank
# sigmamultk=polyval(psigmamultk,tGa);  
  sigmamultk=outer(tGa, 3:0, "^") %*% psigmamultk
  sigma[i]=sigmamultk*mu[i]
  snmean[i]=mu[i] 
  if(mvopt != 0)
   {
# ===== depending on mvopt, delete entries in output vectors =====
    if(mvopt==1) idel= c(6, 7, 10)                             # delete s_i entries except for s_3 and s_4
    if(mvopt==2) idel= 6:10                                    # delete all s_i
    mu= mu[-idel]
    sigma= sigma[-idel]
    alpha= alpha[-idel]
    a1= a1[-idel]
    mu1= mu1[-idel]
    sigma1= sigma1[-idel]
    snmean=snmean[-idel]
   }
  return(list(mu=mu,sigma=sigma,alpha=alpha,a1=a1,mu1=mu1,sigma1=sigma1,snmean=snmean))
 }

# ---------- FUNCTION complogpriorpdf ----------
# complogpriorpdf is an R translation of Alberto Malinverno's MATLAB function.
# it returns value of log-prior pdf of the model vector mv, which contains the
# fundamental Solar system frequencies g_i and s_i, the precession frequency
# k, and the sedimentation rate u.
# - The freqs. g_i and s_i have a skew-normal pdf with parameters from
#   Table 2 of Hoang et al. 2021
# - k has a normal pdf with a mean from a polynomial fit to
#   Figure 6 of Farhat et al. 2022 and a standard deviation from the
#   uncertainties in the calculator of Waltham 2015
# - u has a uniform distribution between umin and umax.
# Input:
# - mv: vector of nm parameters (g_i, s_i, k, and u)
# - prior: list of parameters defining the priors (excluding u), derived from getpriorpdfastropara 
#    (mu,sigma,alpha,a1,mu1,sigma1,snmean)
# - mvopt: defines the content of parameter vector mv (optional, default=0)
#   0: mv = [g_1, ..., g_5, s_1, ..., s_4, s_6, k, u]
#   1: mv = [g_1, ..., g_5, s_3, s_4, k, u]
#   2: mv = [g_1, ..., g_5, k, u]
# - tGa: age in Ga
# - umin, umax: prior bounds on sedimentation rate u
# Output:
# - logpriorpdf, prior pdf for parameter values in mv
# NOTE: if alpha is not zero, the mean of the pdf is
#   delta=alpha/sqrt(1+alpha^2);
#   pdfmean=mu+sigma*delta*sqrt(2/pi);
# NOTE: The secondary mode in the pdfs of g_4 and s_3 is not taken into
# account by the formula above for pdfmean.
#
# Modified to allow for a choice of secular frequencies g_i
# and s_i, as indicated by mvopt. Assumptions:
# - The first 5 elements of mv are always g_1 to g_5
# - The last 2 elements of mv are always k and u, respectively
.complogpriorpdf <- function(mv,prior,mvopt,tGa,umin,umax,iu,ik)
 {
  mu=prior$mu
  sigma=prior$sigma
  alpha=prior$alpha
  a1=prior$a1
  mu1=prior$mu1
  sigma1=prior$sigma1
# make sure the number of model parameters in mv is correct
  nmv = length(mv)
  if (nmv != (length(mu)+1)) stop('Mismatch in the number of parameters in mv')
# compute log-prior pdf
  logpriorpdf = 0
# log-prior pdfs of g_i's and s_i's
  for (i in 1:(ik-1))
   {
    mvipdf <- .skewnormalpdf(mu[i],sigma[i],alpha[i],mv[i])
    if (a1[i] != 0)                                            # include secondary peak
     {
      secpdf <- .skewnormalpdf(mu1[i],sigma1[i],0,mv[i])       # normal pdf, alpha=0
      mvipdf = (1-a1[i])*mvipdf+a1[i]*secpdf
     }
    logpriorpdf = logpriorpdf+log(mvipdf)
   }
# log-prior pdf of precession frequency
  mvipdf <- .skewnormalpdf(mu[ik],sigma[ik],0,mv[ik])          # normal pdf, alpha=0
  logpriorpdf = logpriorpdf+log(mvipdf)
# log-prior pdf of sedimentation rate
  u = mv[iu]
  if ((u>=umin) && (u<=umax))
   {
    logpriorpdf = logpriorpdf+log(1/(umax-umin))
   } else {
    logpriorpdf = -Inf                                         # prior pdf is zero
   }
  return(logpriorpdf)
 }

# ---------- FUNCTION skewnormalpdf ----------
# skewnormalpdf is an R-translation of Alberto Malinverno's MATLAB function.
# compute a skew-normal pdf for parameters mu, sigma, alpha at the value(s) in x
.skewnormalpdf <- function(mu,sigma,alpha,x)
{
  xnorm=((x-mu)/sigma)
  smallphi=(1/sqrt(2*pi))*exp(-0.5*xnorm^2)
# set up error function
  erf <- function(y) 2 * pnorm(y * sqrt(2)) - 1
  bigphi <- 0.5*(1+erf((alpha/sqrt(2))*xnorm))
  snpdf=(2/sigma)*smallphi*bigphi
  return(snpdf)
}

# ---------- FUNCTION arpestim ----------
# arpestim is an R translation of Alberto Malinverno's MATLAB function.
# estimate coefficients, partial ACF, driving noise, and driving noise
# variance for an AR(P) process
# dv is a vector
.arpestim <- function(dv,pdfpara)
 {
 if(pdfpara['ar2'] == 1) {dmean=FALSE ; varmeth=2}
 if(pdfpara['ar2'] == 2) {dmean=TRUE ; varmeth=2}
 if(pdfpara['ar2'] == 3) {dmean=FALSE ; varmeth=1}
 if(pdfpara['ar2'] == 4) {dmean=TRUE ; varmeth=1}
# estimate of AR(P) coefficients (Burg method)
# MATLAB arburg doesn't remove mean
  arp <- ar.burg(dv, aic = FALSE, order.max = pdfpara['P'], demean = dmean, var.method = varmeth)
# residual variance
  varw=arp$var.pred
  phiv=arp$ar                          # AR(P) coefficients
  pacfv=as.vector(arp$partialacf)      # partial ACF
# estimate white noise wv and its ACF
# white noise vector (dv - dvpred), prediction is no good for first P points
  wv=arp$resid
  return(list(phiv=phiv,pacfv=pacfv,wv=wv,varw=varw))
 }

# ---------- FUNCTION comploglikf ----------
# comploglikf is an R translation of Alberto Malinverno's MATLAB function.
# compute log-likelihood of data and precession envelope for TimeOptMCMC.
# This is an "empirical Bayes" version, where the hyperparameters varwv and
# phiv are estimated from the residuals of the spectral fit to the data.
# The log-likelihood for the envelope fit is computed on the basis of the
# residuals of the envelope fit and of the input tauenv.
# Input:
# - mv, vector of model parameters containing
#   - fundamental solar system frequencies g_i (arcsec/yr)
#   - fundamental solar system frequencies s_i (arcsec/yr)
#   - Earth precession frequency k (arcsec/yr)
#   - sedimentation rate u (m/kyr)
# - dat (data frame): 
#   - depth 'vector' (pdfpara.zv of MATLAB)
#   - data 'vector'  (dv of MATLAB), detrended and normalized stratigraphic data
# - fc_dat (data frame)
#   - 'vector' of spatial frequencies
#   - 'vector' of real fourier coeff of dat (pdfpara.ftdv of MATLAB)
#   - 'vector' of imaginary fourier coeff of dat (pdfpara.ftdv of MATLAB)
# - pdfpara contains:
#   - mvopt, which defines the content of mv (0, 1, 2, ...)
#   - iu, the index to sedimentation rate u=mv(iu)
#   - tauenv, variance multiplier of envelope residuals (dimensionless)
#   - P, order of AR(P) process for data residuals (dimensionless)
#   - roll, roll-off rate in Taner filter (dimensionless)
#   - deltafp, delta frequency of precession for Taner filter boundaries
#     (cycles/kyr)
# Output:
# - loglikfdv, log-likelihood of the white noise that drives the assumed
#   AR(P) process in the data fit residuals dv-dvpred
# - loglikfdenv, log-likelihood of the envelope fit residuals assumed
#   uncorrelated with a variance multiplier taudenv
# - dvpred, predicted data
# - denv, envelope of filtered climatic precession in the data
# - denvpred, predicted precession envelope
# - wvest, estimated white noise driving the AR(P) process of data 
#   residuals dv-dvpred 
# - phivest, estimated vector of P coefficients of AR(P) process
#
# Modified to allow for a choice of secular frequencies g_i
# and s_i, as indicated by mvopt. Assumptions:
# - The first 5 elements of mv are always g_1 to g_5
# - The last 2 elements of mv are always k and u, respectively
# - The first 10 columns of Gdv contain sine and cosine terms for 5
#   eccentricity freqs. that are differences g_i - g_j
.comploglikf <- function(mv,dat,fc_dat,pdfpara)
 {
# misc. parameters
  log2pi=log(2*pi)                     # constant in log-likelihood
  N=length(dat[,1])
# compute predicted data, data envelope, predicted envelope, and residuals
  pred <- .compdpred(dat,fc_dat,mv,pdfpara['mvopt'],pdfpara['iu'],pdfpara['roll'],pdfpara['deltafp'])
  dvresid=dat[,2]-pred$dvpred
  denvresid=pred$denv-pred$denvpred
# fit AR(P) model to data residuals, estimate white noise vector wvest
  P=as.numeric(pdfpara['P'])
  arpfit <- .arpestim(dvresid,pdfpara)
# ignore first P points of wv, where prediction is no good  
  wvest=arpfit$wv[(P+1):N]
  varwvest=arpfit$varw
  phivest=arpfit$phiv
  Nw=N-P
# empirical Bayes: log-likelihood of data = log-likelihood of wv  
  loglikfdv=-(Nw/2)*(1+log2pi+log(varwvest))
# empirical Bayes: log-likelihood of precession envelope for a given tauenv
  Neff=N/pdfpara['taudenv']
  vardenvresid=sum(denvresid^2)/N
  loglikfdenv=-(Neff/2)*(1+log2pi+log(vardenvresid))
  return(list(loglikfdv=loglikfdv,loglikfdenv=loglikfdenv,dvpred=pred$dvpred,
              denv=pred$denv,denvpred=pred$denvpred,wvest=wvest,phivest=phivest))
 }

# ---------- FUNCTION compdpred ----------
# compdpred is an R translation of Alberto Malinverno's MATLAB function.
# compute predicted stratigraphic data, precession envelope and predicted
# precession envelope data
# Input:
# - dat (data frame): 
#   - depth 'vector' (pdfpara.zv of MATLAB)
#   - data 'vector'  (dv of MATLAB), detrended and normalized stratigraphic data
# - fc_dat (data frame)
#   - 'vector' of spatial frequencies
#   - 'vector' of real fourier coeff of dat (pdfpara.ftdv of MATLAB)
#   - 'vector' of imaginary fourier coeff of dat (pdfpara.ftdv of MATLAB)
# - mv, vector of model parameters containing
#   - fundamental solar system frequencies g_i (arcsec/yr)
#   - fundamental solar system frequencies s_i (arcsec/yr)
#   - Earth precession frequency k (arcsec/yr)
#   - sedimentation rate u (m/kyr)
# - mvopt, which defines the content of mv (0, 1, 2, ...)
# - iu, the index to sedimentation rate u=mv(iu)
# - tauenv, variance multiplier of envelope residuals (dimensionless)
# - roll, roll-off rate in Taner filter (dimensionless)
# - deltafp, delta frequency of precession for Taner filter boundaries
#   (cycles/kyr)
# Output:
# - dvpred, predicted stratigraphic data
# - denv, precession-band envelope (normalized)
# - denvpred, predicted precession-band envelope (normalized)
# - denvorig, original precession-band envelope (not normalized)
# - bpdv, precession-band filtered data
# - Gdv, matrix G for mv
# - fitcoeffdv, fit coefficients that give dvpred=Gdv*fitcoeffdv
#
# Modified to allow for a choice of secular frequencies g_i
# and s_i, as indicated by mvopt. Assumptions:
# - The first 5 elements of mv are always g_1 to g_5
# - The last 2 elements of mv are always k and u, respectively
# - The first 10 columns of Gdv contain sine and cosine terms for 5
#   eccentricity freqs. that are differences g_i - g_j
.compdpred <- function (dat,fc_dat,mv,mvopt,iu,roll,deltafp)
 {
  zv=dat[,1]
  dv=dat[,2]
# misc. parameters
  fconv=1000/(360*60*60)               # conversion from arcsec/yr to cycles/kyr
  Ne=5*2                               # number of columns in Gdv with eccentricity frequencies
  Nd=length(dv)
# compute ages (tv) for data series, with 0 at the top
  u=mv[iu]
  tv=(zv-zv[1])/u                      # kyr
#  deltat=tv[2]-tv[1]
# set up bandpass filter to compute data envelope from the Hilbert transform
# climatic precession frequencies for this iteration
  pv <- .compeopfreq(mv,mvopt)$pv
  pv=fconv*pv
  bpfmin=min(pv)-deltafp
  bpfmax=max(pv)+deltafp
#  convert frequencies in fc_dat from spatial to temporal, given present sedimentation rate
  fc_dat2=fc_dat
  fc_dat2[1]=fc_dat2[1]*u
# apply taner filter to fourier coeffs
  bpdat <- tanerFC(fc_dat2,npts=Nd,flow=bpfmin,fhigh=bpfmax,roll=roll,output=1,genplot=F,verbose=F)
  bpdv = bpdat[,2]                     # return values
# envelope of climatic precession in the data
#  note: original Matlab hilbert implementation (M&M 2024) is better reproduced with padfac=1, but padfac=2
#  is used here to address edge effects
  denvorig <- hilbert(bpdat,padfac=2,demean=T,detrend=F,addmean=F,genplot=F,check=F,verbose=F)[,2]
# standardize and save envelope to vector (dv is already standardized)
  denv=denvorig
  denv=denv-mean(denv)
  denv=denv/sd(denv)
# compute data matrix Gdv and envelope matrix Gdenv
  Gmats <- .compGmats(mv,mvopt,tv)
  Gdv=Gmats$Gdv
  Gdenv=Gmats$Gdenv
# compute dvpred=Gdv*mdv
# GdvTGdv=Gdv'*Gdv;
  GdvTGdv=t(Gdv) %*% Gdv
# fitcoeffdv=GdvTGdv\(Gdv'*dv); % least-squares fit coeffs.
  fitcoeffdv <- solve(GdvTGdv,(t(Gdv) %*% dv))                 # least-squares fit coeffs.
  dvpred=Gdv %*% fitcoeffdv            # predicted data
# compute denvpred=Gdenv*mdenv
  GdenvTGdenv=GdvTGdv[1:Ne,1:Ne]       # no need to recalculate GdenvTGdenv
# mdenv=GdenvTGdenv\(Gdenv'*denv); % least-squares fit coeffs.
  mdenv <- solve(GdenvTGdenv,(t(Gdenv) %*% denv))              #  least-squares fit coeffs.
  denvpred=Gdenv %*% mdenv             # predicted envelope
  return(list(dvpred=dvpred,denv=denv,denvpred=denvpred,denvorig=denvorig,bpdv=bpdv,Gdv=Gdv,fitcoeffdv=fitcoeffdv))
 }

# ---------- FUNCTION compGmats ----------
# compGmats is an R translation of Alberto Malinverno's MATLAB function.
# compute G matrices for spectral and envelope fit in TimeOpt
# Input:
# - mv, vector of 12 model parameters containing
#   - five fundamental solar system frequencies g_i (arcsec/yr)
#   - five fundamental solar system frequencies s_i (arcsec/yr)
#   - Earth precession frequency k (arcsec/yr)
#   - sedimentation rate u (m/kyr)
# - tv, vector of ages (kyr)
# Output:
# - Gdv, matrix G for spectral fit
# - Gdenv, matrix G for envelope fit
#
# Modified to allow for a choice of secular frequencies g_i
# and s_i, as indicated by mvopt. Assumptions:
# - The first 5 elements of mv are always g_1 to g_5
# - The last 2 elements of mv are always k and u, respectively
# - The first 10 columns of Gdv contain sine and cosine terms for 5
#   eccentricity freqs. that are differences g_i - g_j
.compGmats <- function (mv,mvopt,tv)
 {
# parameters
  fconv=1000/(360*60*60)               # conversion from arcsec/yr to cycles/kyr
  Ng=5                                 # number of fundamental solar system frequencies
  Ns=5
  Nd=length(tv)                        # number of data points
# eccentricity, obliquity, and climatic precession frequencies
  eop <- .compeopfreq(mv,mvopt)
  ev=fconv*eop$ev                      # convert arcsec/yr to cycles/kyr
  ov=fconv*eop$ov                      # ov can be NA
  pv=fconv*eop$pv
# allocate spectral fit matrix Gdv in dvpred=Gdv*mdv
  Ne=length(ev)
  No=length(na.omit(data.frame(ov))[,1])                       # may be zero (if no obliquity used, ov=NA returned from .compeopfreq)
  Np=length(na.omit(data.frame(pv))[,1])          
  Ncol=2*(Ne+No+Np)                    # total number of columns in Gdv
  Gdv=double(Nd*Ncol)
  Gdv=matrix(0, Nd, Ncol)
# fill matrix Gdv
  j=1                                  # index of column in Gdv
  for (ie in 1:Ne)                     # Ne eccentricity freqs. 
   {
    arg=2*pi*ev[ie]*tv
    Gdv[,j] <- cos(arg)
    Gdv[,j+1] <- sin(arg)
    j=j+2
   }
  if(No > 0)                           # No obliquity freqs. (No may be zero)
   {
    for (io in 1:No)                     
     {
      arg=2*pi*ov[io]*tv
      Gdv[,j]=cos(arg)
      Gdv[,j+1]=sin(arg)
      j=j+2
     }
   }  
  if(Np > 0)                           # Np climatic prec. freqs. (Np may be zero)
   {
    for (ip in 1:Np)                     
     {
      arg=2*pi*pv[ip]*tv
      Gdv[,j]=cos(arg)
      Gdv[,j+1]=sin(arg)
      j=j+2
     }
   }  
# envelope fit matrix Gdenv in denvpred=Gdenv*mdenv
  Gdenv=Gdv[,1:(2*Ne)]                 # Gdenv only has Ne eccentricity freqs.
  return(list(Gdv=Gdv,Gdenv=Gdenv))
 }

# ---------- FUNCTION comppriorpdf1 ----------
# comppriorpdf1 is an R translation of Alberto Malinverno's MATLAB function.
# return value of prior pdf of a single model parameter (one of the
# fundamental Solar system frequencies g_i or s_i, the precession
# frequency k,or the sedimentation rate u). 
# - The freqs. g_i and s_i have a skew-normal pdf with parameters from
#   Table 2 of Hoang et al. 2021
# - k has a normal pdf with a mean from a polynomial fit to 
#   Figure 6 of Farhat et al. 2022 and a standard deviation from the
#   uncertainties in the calculator of Waltham 2015
# - u has a uniform distribution between umin and umax
# Input:
# - mvstr: text string of model parameter (e.g., g_3, s_4, k)
# - mvval: the prior pdf is calculated for the values in vector mvval 
# - mvopt, which defines the content of parameter vector mv
# - tGa: age in Ga
# Output:
# - priorpdf, prior pdf for parameter value(s) in mvval
#
# Modified to allow for a choice of secular frequencies g_i
# and s_i, as indicated by mvopt. Assumptions:
# - The first 5 elements of mv are always g_1 to g_5
# - The last 2 elements of mv are always k and u, respectively
# - The first 10 columns of Gdv contain sine and cosine terms for 5
#   eccentricity freqs. that are differences g_i - g_j
.comppriorpdf1 <- function(mvstr,mvval,mvopt,tGa,umin,umax)
 {
# compute prior pdf of parameter in mvstr
  if (mvstr == 'u')                    # prior pdf of sed. rate
   {
    priorpdf=double(length(mvval))
    priorpdf[mvval >= umin & mvval <= umax] = 1 / (umax-umin)
   }else{                              # prior pdf of one of the astronomical frequencies
    if(tGa>0)
     {    
      imv <- .getmvindex(mvstr,mvopt)
      prior <- .getpriorpdfastropara(tGa,mvopt)
      mu = prior$mu
      sigma = prior$sigma
      alpha = prior$alpha
      a1 = prior$a1
      mu1 = prior$mu1
      sigma1 = prior$sigma1                    
      priorpdf <- .skewnormalpdf(mu[imv],sigma[imv],alpha[imv],mvval)
      if (a1[imv] != 0)                # include secondary peak
       {
        secpdf <- .skewnormalpdf(mu1[imv],sigma1[imv],0,mvval) # normal pdf, alpha=0
        priorpdf=(1-a1[imv])*priorpdf+a1[imv]*secpdf
       }
     }else{
        priorpdf=double(length(mvval))
     }  
   }
  return(priorpdf)
 }

# ---------- FUNCTION comppdfpara ----------
# comppdfpara is an R translation of Alberto Malinverno's MATLAB function.
# compute parameters of an unnormalized input pdf
# Input:
# - x: vector of evenly spaced x-values
# - pdfux: unnormalized pdf of x 
# - alpha: value of confidence/credible interval (0-1)
# Output:
# - xmean: mean of x
# - xsdev: standard deviation of x
# - xlow, xhigh: boundaries of alpha confidence/credible interval
# - pdfx: normalized pdf of x 
# - cdfx: cumulative distribution function of x 
.comppdfpara <- function(x,pdfux,alpha)
 {
# compute normalized pdf and cdf
  pdfx=pdfux/sum(pdfux)                                        # pdfx sums to 1
  cdfx=cumsum(pdfx)                                            # cdfx is between 0 and 1
# compute mean and standard deviation
  xmean=sum(x*pdfx)
  xvar=sum(((x-xmean)^2)*pdfx)
  xsdev=sqrt(xvar)
# compute alpha confidence/credible interval
  cdflow=(1-alpha)/2
  cdfhigh=1-cdflow
#  ilow=find(cdfx>cdflow,1,'first');
  ilow=which(cdfx>cdflow)[1]
  if (ilow>1) 
   {
    xlow <- approx(c(cdfx[ilow-1],cdfx[ilow]), c(x[ilow-1], x[ilow]), cdflow, method = "linear")$y
   }else{
    xlow=x[1]
    cat("comppdfpara: WARNING: pdfx did not contain low",100*alpha, "% bound\n")
   }
  ihigh=tail(which(cdfx<cdfhigh),n=1)
  if (ihigh<length(cdfx))
   {
    xhigh <- approx(c(cdfx[ihigh],cdfx[ihigh+1]), c(x[ihigh], x[ihigh+1]), cdfhigh, method = "linear")$y
   }else{
    xhigh=x[length(x)]
    cat("comppdfpara: WARNING: pdfx did not contain high", 100*alpha,"% bound\n")
   }
# correct pdf so it integrates to 1 (rather than sum to 1)
  deltax=x[2]-x[1]                                             # x must be evenly spaced
  pdfx=pdfx/deltax
  return(list(xmean=xmean,xsdev=xsdev,xlow=xlow,xhigh=xhigh,pdfx=pdfx,cdf=cdfx))
 }
 
# ---------- FUNCTION comppvalue ----------
# comppvalue is an R translation of Alberto Malinverno's MATLAB function. 
# returns the p-value (0 to 1) given
# - xdata = value of a data statistic (e.g., R^2)
# - xsim = Nsim values of the same statistic computed for random simulated
#   data.
# The p-value is the frequency of xsim values >= xdata. Also returns the
# number of xsim values >= xdata as pvalcount.
.comppvalue <- function(xdata,xsim)
 {
  Nsim=length(xsim)
  pvalindex=which(xsim>=xdata)
  pvalcount=length(pvalindex)
  if (pvalcount>0)
   { 
    pvalue=pvalcount/Nsim
   }else{ 
    pvalue=1/Nsim
   }   
   return(list(pvalue,pvalcount))
  } 

# ---------- FUNCTION comprsq ----------
# comprsq is an R translation of Alberto Malinverno's MATLAB function.
# returns the R^2 value for the given model parameters and data,
# written to minimize calculations when it is called in the Monte Carlo 
# simulation to assess significance
# Input:
# - dv: observed data vector
# - sumsqrdv: sum of squares of vector dv (so it is not repeatedly computed
#   in Monte Carlo simulations)
# - G: matrix that gives predicted data as in dvpred=G*afitv
# - afitv: least-squares coeffs.
# - jafitv: if present, zero out entries in afitv that are not in the
#   indices in jaftiv; this allows for computing R^2 when only some cycles
#   are included (eccentricity, obliquity, or clim.prec.)
# Output:
# - rsq: R^2 between dvpred=G*afitv and dv
.comprsq <- function(dv,sumsqrdv,G,afitv,jafitv=NULL)
 {
  bfitv=afitv
  if (!is.null(jafitv)) 
   {
     bfitv[setdiff(1:length(bfitv),jafitv)] = 0        # zero out elements of bfitv whose indices are not in jafitv
   }
  dvpred= G %*% bfitv
  dvresid= dv - dvpred
  rsq=(sumsqrdv-sum(dvresid^2))/sumsqrdv
  return(rsq)
 }
 
# ---------- FUNCTION plotmvdatafit ----------
# plotmvdatafit is an R translation of Alberto Malinverno's MATLAB function.
# function that plots the fit to the data for a given model parameter
# vector mv in the current figure.
# Input:
#  dat (data frame) 
#    - depth 'vector' (zv of MATLAB)
#    - data 'vector' (dv of MATLAB), detrended and normalized
#  fc_dat (data frame)
#    - 'vector' of spatial frequencies
#    - 'vector' of real fourier coeff of dat (pdfpara.ftdv of MATLAB)
#    - 'vector' of imaginary fourier coeff of dat (pdfpara.ftdv of MATLAB)
#  mv (model vector) <---- Typically we want to use the MAP values
#    - The first 5 elements of mv are always g_1 to g_5
#    - The last 2 elements of mv are always k and u, respectively
#    - The first 10 columns of Gdv will contain sine and cosine terms for 5
#      eccentricity freqs. that are differences g_i - g_j
# - mvopt: defines the content of parameter vector mv
# - nolikfdenv: the likelihood for the envelope fit is not used (optional,
#   default=0)
# Output:
# - rsqdatafit: signal/noise ratio of data fit
# - rsqdenvfit: signal/noise ratio of precession envelope fit
#
# Modified to allow for a choice of secular frequencies g_i
# and s_i, as indicated by mvopt. Assumptions:
# - The first 5 elements of mv are always g_1 to g_5
# - The last 2 elements of mv are always k and u, respectively
# - The first 10 columns of Gdv contain sine and cosine terms for 5
#   eccentricity freqs. that are differences g_i - g_j
.plotmvdatafit <- function(dat,fc_dat,mv,mvopt,iu,roll,deltafp,nolikfdenv=0)
 {
  zv=dat[,1]
  dv=dat[,2]
  Ndv=length(dv)
# parameters
  evcolor="purple"
  ovcolor="orange"
  pvcolor="seagreen"

# compute predicted data, predicted envelope, and R^2 values
  pred <- .compdpred(dat,fc_dat,mv,mvopt,iu,roll,deltafp)
  dvpred = pred$dvpred
  denv = pred$denv
  denvpred = pred$denvpred
  denvorig = pred$denvorig
  bpdv = pred$bpdv
  sumsqrdv=sum(dv^2)                   # assume zero mean in standardized dv
  sumsqrdvresid=sum((dv-dvpred)^2)
  rsqdatafit=(sumsqrdv-sumsqrdvresid)/sumsqrdv
  sumsqrdenv=sum((denv-mean(denv))^2)
  sumsqrdenvresid=sum((denv-denvpred)^2)
  rsqdenvfit=(sumsqrdenv-sumsqrdenvresid)/sumsqrdenv           # assume zero mean

# compute original (unnormalized) predicted envelope for plotting
  mu=mean(denvorig)
  sigma=sd(denvorig)
  denvpredorig=mu+denvpred*sigma

# compute periodogram power
  tv=(zv-zv[1])/mv[iu]                 # time vector in kyr
  datT=data.frame(cbind(tv,dv))
# increasing padding from 2 to 8  
  pgdat=periodogram(datT,padfac=8,demean=T,detrend=F,output=1,nrm=1,genplot=F,check=F,verbose=F)
  fv=pgdat[,1]
  pgdv=pgdat[,3]

# compute predicted eccentricity, obliquity, climatic precession freqs. for
# the MAP model parameter vector
  eop <- .compeopfreq(mv,mvopt,label=TRUE)
  ev=eop$ev
  ov=eop$ov
  pv=eop$pv
  evlabel=eop$evlabel
  ovlabel=eop$ovlabel
  pvlabel=eop$pvlabel
#  if (any(is.na(pv))) stop('**** ERROR: Need a vector of climatic precession frequencies. TERMINATING NOW!')
  fconv=1000/(360*60*60)               # conversion from arcsec/yr to cycles/kyr
  ev=fconv*ev
  ov=fconv*ov
  pv=fconv*pv
  maxfplot=1.25*max(pv)                # max. frequency to plot in periodogram
  maxpgdv=1.2*max(pgdv[fv<maxfplot])   # max. power to plot in periodogram, in given frequency range
  minpgdv=min(pgdv[fv<maxfplot])       # min. power to plot in periodogram, in given frequency range
  midpgdv=(minpgdv+maxpgdv)/2
  
# compute frequency response of Taner filter used to bandpass filter
# climatic precession
  bpfmin=min(pv)-deltafp
  bpfmax=max(pv)+deltafp
# call taner filter again, due to new padding (for comparison with pgdat)
  bpfresp <- taner(datT,padfac=8,flow=bpfmin,fhigh=bpfmax,roll=roll,output=2,genplot=F,verbose=F)
# mid-value of depth for plotting R^2 value
  zmid=(zv[1]+zv[Ndv])/2

# ===== plot fit to the data =====
  dev.new(title=c("MAP data fit"),height = 7.3, width = 5.9, units = "in")
  layoutA <- layout(matrix(c(1,2,3), 3, 1)) 
# margins: bottom, left, top, right
  par(mar = c(3.5, 3.6, 2.1, 1))
#         title, text label, axis_line
  par(mgp = c(3, 0.5, 0))
# plot data and predicted data
  plot(zv,dv,type="l",xlab="",ylab="",main="MAP data fit",ylim=c(min(dv,dvpred),max(dv,dvpred)),bty="n")
  axis(side=1, at=c(zv[1],zv[Ndv]), labels=c("",""), lwd.ticks=0)
  lines(zv,dvpred,col="red",lwd=2)
  mtext(c("Depth (m)"),side=1,line=1.6,cex=0.8)
  mtext(c("Standardized data"),side=2,line=1.6,cex=0.8)
  legend(x="bottomright",legend=c('Data','Predicted data'),col=c("black","red"),lty=c(1,1),lwd=c(1,2),bg="white",cex=1)
  mtext(bquote(r^2==.(round(rsqdatafit,digits=2))),side=1,line=-1,cex=0.9)

# plot envelope and predicted data 
# margins: bottom, left, top, right
  par(mar = c(3.2, 3.6, 3.1, 1))
  plot(zv,bpdv,type="l",xlab="",ylab="",main="",col="gray",ylim=c(min(bpdv,denvorig,denvpredorig),max(bpdv,denvorig,denvpredorig)),bty="n")
  axis(side=1, at=c(zv[1],zv[Ndv]), labels=c("",""), lwd.ticks=0)
  lines(zv,denvorig,col="black",lwd=2)
  lines(zv,denvpredorig,col="red",lwd=2)
  mtext(c("Depth (m)"),side=1,line=1.6,cex=0.8)
  mtext(c("Standardized data"),side=2,line=1.6,cex=0.8)
  legend(x="bottomright",legend=c('Bandpass filtered data','Envelope','Predicted envelope'),col=c("gray","black","red"),lty=c(1,1,1),lwd=c(1,2,2),bg="white",cex=1)
# cross out envelope result if not used
  if (nolikfdenv == 1) 
   {
    lines(c(zv[1],zv[Ndv]),c(min(bpdv,denvorig,denvpredorig),max(bpdv,denvorig,denvpredorig)),col="#00000046",lwd=10)
    lines(c(zv[1],zv[Ndv]),c(max(bpdv,denvorig,denvpredorig),min(bpdv,denvorig,denvpredorig)),col="#00000046",lwd=10)
   }
# plot second (time) axis on top
  par(new = TRUE)
  plot(tv,dv,type="n",xlab="",ylab="",main="",xaxt="n",yaxt="n",bty="n")
  axis(side=3,line=0,col="royalblue",cex=0.7,col.axis="royalblue",font=3)
  axis(side=3,line=0,col="royalblue",at=c(tv[1],tv[Ndv]), labels=c("",""), lwd.ticks=0)
  mtext(c("Time (kyr)"),side=3,line=1.5,cex=0.7,col="royalblue",font=3)
  mtext(bquote(r^2==.(round(rsqdenvfit,digits=2))),side=1,line=-1,cex=0.9)

# margins: bottom, left, top, right
  par(mar = c(3.1, 3.6, 1.1, 1))
# plot periodogram, filter response, and predicted frequencies
# set-up plot
  plot(0,0,cex=0,type="n",xlab="",ylab="",main="",xlim=c(0,maxfplot),ylim=c(minpgdv,1.05*maxpgdv))
  mtext(c("Frequency (cycles/kyr)"),side=1,line=1.6,cex=0.8)
  mtext(c("Data periodogram"),side=2,line=1.6,cex=0.8)
# plot filter
  polygon(bpfresp[,1],bpfresp[,2]*maxpgdv,col="#FFFF005A",border=NA)
# plot and label eccentricity, obliquity, climatic precession freqs.
  Ne=length(ev)
  No=length(na.omit(data.frame(ov))[,1])                       # may be zero (if no obliquity used, ov=NA returned from .compeopfreq)
  Np=length(pv)
  for (i in seq(Ne,1,-1))
   {
    if(i%%2 != 0) ynow=mean(c(midpgdv,maxpgdv))
    if(i%%2 == 0) ynow=midpgdv
    ycord=c(minpgdv,ynow)
    lines(x=c(ev[i],ev[i]),y=ycord,col=evcolor,lty=3,lwd=1.5)
    text(x=ev[i],y=ynow,evlabel[i],srt=90,font=3,col=evcolor,pos=3,offset=1.25)
   }
  text((max(ev)+min(ev))/2,1.1*maxpgdv,"Eccentricity",col=evcolor,cex=1.5,pos=1)

  if (No>0)
   {
    for (i in seq(No,1,-1))
     {
      if(i%%2 != 0) ynow=mean(c(midpgdv,maxpgdv))
      if(i%%2 == 0) ynow=midpgdv
      ycord=c(minpgdv,ynow)
      lines(x=c(ov[i],ov[i]),y=ycord,col=ovcolor,lty=3,lwd=1.5)
      text(x=ov[i],y=ynow,ovlabel[i],srt=90,font=3,col=ovcolor,pos=3,offset=1.25)
     }
    text((max(ov)+min(ov))/2,1.1*maxpgdv,"Obliquity",col=ovcolor,cex=1.5,pos=1)
   }

  if (Np>0)
   {
    for (i in seq(Np,1,-1))
     {
      if(i%%2 != 0) ynow=mean(c(midpgdv,maxpgdv))
      if(i%%2 == 0) ynow=midpgdv
      ycord=c(minpgdv,ynow)
      lines(x=c(pv[i],pv[i]),y=ycord,col=pvcolor,lty=3,lwd=1.5)
      text(x=pv[i],y=ynow,pvlabel[i],srt=90,font=3,col=pvcolor,pos=3,offset=1.25)
     }
    text((max(pv)+min(pv))/2,1.1*maxpgdv,"Precession",col=pvcolor,cex=1.5,pos=1)
   } 

# plot periodogram
  lines(fv,pgdv)
  return(invisible(list(rsqdatafit=rsqdatafit,rsqdenvfit=rsqdenvfit)))
 } 
#######################################################################################
