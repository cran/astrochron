### This function is a component of astrochron: An R Package for Astrochronology
### Copyright (C) 2026 Stephen R. Meyers
###
##################################################################################
### function timeOptBMCMCplot - (SRM: Feb. 21, 2026)
##################################################################################

### NOTE: if you aren't calling this function from within timeOptBMCMC, you will need
#         to explicitly detrend if that option was selected in timeOptBMCMC, 
#         and standardize the data too. Or use the 'dat' returned from timeOptBMCMC,
#         which has already been prepared.

timeOptBMCMCplot <- function(dat,mcmcres,pdfpara,nburnin=NULL,fc_dat=NULL) 
 {
  if(pdfpara['prioronly'] == 1) test=TRUE
  if(pdfpara['prioronly'] == 0) test=FALSE
# mcmcres$mv has log pdf in first column, which results in
# adding or subtracting 1 to ik, iu
  ik=pdfpara['ik'] + 1                 # mv values start in column 2 of mcmcres$mv
  iu=pdfpara['iu'] + 1                 # mv values start in column 2 of mcmcres$mv
  mvopt=pdfpara['mvopt']

  niter=length(mcmcres$mv[,1])
  nbatch=length(mcmcres$pacc[,1])
  nmv=dim(mcmcres$mv)[2]-1
  mvlabels=.getmvlabels(mvopt)

  alpha=0.95                           # P of credible interval

  pngwidth= 800
  pngheight= 800
  pngres= 130
  cexset= 0.75

  if(is.null(fc_dat)) fc_dat <- periodogram(dat,padfac=2,demean=T,detrend=F,output=2,nrm=0,genplot=F,check=F,verbose=F)

  if(is.null(nburnin))
   {
    if(!pdfpara['prioronly'] == 1)
     {
# burn in detection, using median value from second half of sample pdf
      logpdfburnin=median(mcmcres$mv[(niter/2):niter,1])
# identify the first value to go above it, use that as the start of the burn-in  
      nburnin=which(mcmcres$mv[,1]>logpdfburnin)[1]
      cat(" * MCMC chain burn-in: discard",nburnin-1, "samples.\n")
     }else{
      nburnin=1
     }
   }

# determine MAP sampled mv, accounting for burn-in
  itermax = which.max(mcmcres$mv[nburnin:niter,1]) + nburnin - 1
# note that in mcmcres$mv, the mvs start with second column  
  mvmap=mcmcres$mv[itermax,]

# calculate posterior statistics of u and k, accounting for burn-in
  kmap=mvmap[ik]
  umap=100*mvmap[iu]                   # cm/kyr
  ksample=mcmcres$mv[nburnin:niter,ik]
  usample=100*mcmcres$mv[nburnin:niter,iu]                     # cm/kyr
  kpostmean=mean(ksample)
  kpostsdev=sd(ksample)
  upostmean=mean(usample)
  upostsdev=sd(usample)
# bounds of specified credible interval
  ilow=round(length(ksample)*(1-alpha)/2)
  ihigh=length(ksample)-ilow
# also define 99.9% credible interval, for plotting constraint on g4-3, s4-s3
  ilow999=round(length(ksample)*(1-0.999)/2)
  ihigh999=length(ksample)-ilow999
  ksamplesort=sort(ksample)
  kpostlow=ksamplesort[ilow]
  kposthigh=ksamplesort[ihigh]
  usamplesort=sort(usample)
  upostlow=usamplesort[ilow]
  uposthigh=usamplesort[ihigh]

# calculate periods of g4 - g3 and bounds of credible interval, accounting
# for burnin
  pconv=(360*60*60)*1e-6               # conversion from period of yr/arcsec to Myr/cycle
  ig3 <- .getmvindex('g3',mvopt) + 1                 # mv values start in column 2 of mcmcres$mv
  ig4 <- .getmvindex('g4',mvopt) + 1                 # mv values start in column 2 of mcmcres$mv
  g4mg3=pconv*(mcmcres$mv[nburnin:niter,ig4]-mcmcres$mv[nburnin:niter,ig3])^(-1)
  g4mg3samplesort=sort(g4mg3)
  g4mg3postlow=g4mg3samplesort[ilow]
  g4mg3posthigh=g4mg3samplesort[ihigh]

  if (mvopt <2)  
   {
    is3 <- .getmvindex('s3',mvopt) + 1                # mv values start in column 2 of mcmcres$mv
    is4 <- .getmvindex('s4',mvopt) + 1                # mv values start in column 2 of mcmcres$mv
    s4ms3=pconv*(mcmcres$mv[nburnin:niter,is4]-mcmcres$mv[nburnin:niter,is3])^(-1)
    s4ms3samplesort=sort(s4ms3)
    s4ms3postlow=s4ms3samplesort[ilow]
    s4ms3posthigh=s4ms3samplesort[ihigh]
   } 

# report posterior statistics 
   cat("\n========= Marginal posterior PDF =================================\n")
   cat(sprintf('  \tMAP  ,\tmean  ,\tst.dev  ,\t%g%% interval\n',100*alpha))
   cat(sprintf('u:\t%.6f  ,\t%.6f  ,\t%.6f  ,\t%.6f - %.6f\n',umap,upostmean,upostsdev,upostlow,uposthigh))
   cat(sprintf('k:\t%.6f  ,\t%.6f  ,\t%.6f  ,\t%.6f - %.6f\n',kmap,kpostmean,kpostsdev,kpostlow,kposthigh))
   cat(sprintf('g4-g3:\t------- ,\t------- ,\t------ ,\t%.6f - %.6f\n',g4mg3postlow,g4mg3posthigh))
   if(mvopt<2) cat(sprintf('s4-s3:\t------- ,\t------- ,\t------ ,\t%.6f - %.6f\n',s4ms3postlow,s4ms3posthigh))


# PLOT 1: MCMC diagnostic plots
  if(!test)
   {
    dev.new(title=c("Adaptive MCMC diagnostics"))
    if(mvopt == 0 || mvopt == 1) layoutA <- layout(matrix(c(1,1,1,2,3,4,5,6,7), 3, 3, nrow=3)) 
    if(mvopt == 2) layoutA <- layout(matrix(c(1,1,2,3,4,5), 2, 3, nrow=3)) 
#          bottom, left, top, right
    par(mar = c(3.2, 3.6, 2.1, 1))
# logpdf versus iteration
    plot(mcmcres$mv[,1],xlab="",ylab="",main="",cex=cexset,type="l")
    abline(v=nburnin,col="red",lwd=2,lty=3)   
    mtext("Sample number", side=1,line=2.3)
    mtext("Log-posterior PDF", side=2,line=2.2)
    mtext("MCMC sampling",side=3,line=0.5)
#          bottom, left, top, right
    par(mar = c(3.2, 3.6, 1.5, 1))
# sigmamvcand proposal/start versus iteration
    pacc=mcmcres$pacc
    sigmamvR=mcmcres$sigmamvcand
    sigmamvR <- lapply(1:nmv, function(cols) sigmamvR[,cols] / sigmamvR[1,cols] )
    sigmamvR = matrix(unlist(sigmamvR), ncol = nmv, byrow = FALSE)

# plot g_i
    plot((1:nbatch)*50,sigmamvR[,1],xlab="",ylab="",main="",cex=cexset,col="blue",type="l",ylim=c(min(sigmamvR[,1:5]),max(sigmamvR[,1:5])))
    lines((1:nbatch)*50,sigmamvR[,2],col="red")
    lines((1:nbatch)*50,sigmamvR[,3],col="yellow")
    lines((1:nbatch)*50,sigmamvR[,4],col="purple")
    lines((1:nbatch)*50,sigmamvR[,5],col="green")
    mtext("Sample number", side=1,line=2.3)
    mtext("sigma proposal/start", side=2,line=2.2)
    if(max(sigmamvR[nbatch,1:5]) > 1) legloc="bottomright"
    if(max(sigmamvR[nbatch,1:5]) < 1) legloc="topright"
    legend(x=legloc,legend=c('g1','g2','g3','g4','g5'),col=c("blue","red","yellow","purple","green"),lty=c(1,1,1,1,1),lwd=c(1,1,1,1,1),cex=0.7)

# plot s_i
    if(mvopt == 0)                     # include all s_i
     {
      plot((1:nbatch)*50,sigmamvR[,6],xlab="",ylab="",main="",cex=cexset,col="blue",type="l",ylim=c(min(sigmamvR[,6:10]),max(sigmamvR[,6:10])))
      lines((1:nbatch)*50,sigmamvR[,7],col="red")
      lines((1:nbatch)*50,sigmamvR[,8],col="yellow")
      lines((1:nbatch)*50,sigmamvR[,9],col="purple")
      lines((1:nbatch)*50,sigmamvR[,10],col="green")
      if(max(sigmamvR[nbatch,6:10]) > 1) legloc="bottomright"
      if(max(sigmamvR[nbatch,6:10]) < 1) legloc="topright"
      legend(x=legloc,legend=c('s1','s2','s3','s4','s5'),col=c("blue","red","yellow","purple","green"),lty=c(1,1,1,1,1),lwd=c(1,1,1,1,1),cex=0.7)
      mtext("Sample number", side=1,line=2.3)
      mtext("sigma proposal/start", side=2,line=2.2)
     } 
    if(mvopt == 1)                     # include s_3 and s_4
     {
      plot((1:nbatch)*50,sigmamvR[,6],xlab="",ylab="",main="",cex=cexset,col="blue",type="l",ylim=c(min(sigmamvR[,6:7]),max(sigmamvR[,6:7])))
      lines((1:nbatch)*50,sigmamvR[,7],col="red")
      if(max(sigmamvR[nbatch,6:7]) > 1) legloc="bottomright"
      if(max(sigmamvR[nbatch,6:7]) < 1) legloc="topright"
      legend(x=legloc,legend=c('s3','s4'),col=c("blue","red"),lty=c(1,1),lwd=c(1,1),cex=0.7)
      mtext("Sample number", side=1,line=2.3)
      mtext("sigma proposal/start", side=2,line=2.2)
     } 

# if mvopt == 2, skip plot
#    if(mvopt == 2)                     # no s_i
#     {
#      plot((1:nbatch)*50,rep(0,nbatch),type="n",xlab="",ylab="",main="",cex=cexset,col="blue",ylim=c(0,1))
#      mtext("Sample number", side=1,line=2.3)
#      mtext("sigma proposal/start", side=2,line=2.2)
#     } 

# plot k, u (ik and iu were set for mcmcres$mv)
    plot((1:nbatch)*50,sigmamvR[,ik-1],xlab="",ylab="",main="",cex=cexset,col="blue",type="l",ylim=c(min(sigmamvR[,(ik-1):(iu-1)]),max(sigmamvR[,(ik-1):(iu-1)])))
    lines((1:nbatch)*50,sigmamvR[,iu-1],col="red")
    mtext("Sample number", side=1,line=2.3)
    mtext("sigma proposal/start", side=2,line=2.2)
    if(max(sigmamvR[nbatch,(ik-1):(iu-1)]) > 1) legloc="bottomright"
    if(max(sigmamvR[nbatch,(ik-1):(iu-1)]) < 1) legloc="topright"
    legend(x=legloc,legend=c('k','u'),col=c("blue","red"),lty=c(1,1),lwd=c(1,1),cex=0.7)

# plot g_i
    plot((1:nbatch)*50,pacc[,1],xlab="",ylab="",main="",cex=cexset,col="blue",type="l",ylim=c(min(pacc),max(pacc)))
    lines((1:nbatch)*50,pacc[,2],col="red")
    lines((1:nbatch)*50,pacc[,3],col="yellow")
    lines((1:nbatch)*50,pacc[,4],col="purple")
    lines((1:nbatch)*50,pacc[,5],col="green")
    abline(h=0.4,col="gray",lwd=2,lty=3)   
    mtext("Sample number", side=1,line=2.3)
    mtext("pacc", side=2,line=2.2)
# plot s_i
    if(mvopt == 0)                     # include all s_i
     {
      plot((1:nbatch)*50,pacc[,6],xlab="",ylab="",main="",cex=cexset,col="blue",type="l",ylim=c(min(pacc),max(pacc)))
      lines((1:nbatch)*50,pacc[,7],col="red")
      lines((1:nbatch)*50,pacc[,8],col="yellow")
      lines((1:nbatch)*50,pacc[,9],col="purple")
      lines((1:nbatch)*50,pacc[,10],col="green")
      abline(h=0.4,col="gray",lwd=2,lty=3)   
      mtext("Sample number", side=1,line=2.3)
      mtext("pacc", side=2,line=2.2)
     }
    if(mvopt == 1)                     # include s_3 and s_4
     {
      plot((1:nbatch)*50,pacc[,6],xlab="",ylab="",main="",cex=cexset,col="blue",type="l",ylim=c(min(pacc),max(pacc)))
      lines((1:nbatch)*50,pacc[,7],col="red")
      mtext("Sample number", side=1,line=2.3)
      mtext("pacc", side=2,line=2.2)
     } 

# if mvopt == 2, skip plot
#    if(mvopt == 2)                     # no s_i
#     {
#      plot((1:nbatch)*50,rep(0,nbatch),type="n",xlab="",ylab="",main="",cex=cexset,col="blue",ylim=c(0,1))
#      mtext("Sample number", side=1,line=2.3)
#      mtext("pacc", side=2,line=2.2)
#     }

# plot k, u (ik and iu were set for mcmcres$mv)
    plot((1:nbatch)*50,pacc[,ik-1],xlab="",ylab="",main="",cex=cexset,col="blue",type="l",ylim=c(min(pacc),max(pacc)))
    lines((1:nbatch)*50,pacc[,iu-1],col="red")
    abline(h=0.4,col="gray",lwd=2,lty=3)   
    mtext("Sample number", side=1,line=2.3)
    mtext("pacc", side=2,line=2.2)
   }


# PLOT 2: posterior histograms with priors
# 12 plots for g, s, k, u
  dev.new(title=c("Posterior & prior PDFs: Fundamental frequencies, k, u"))
  layoutB <- layout(matrix(c(1,2,3,4,5,6,7,8,9,10,11,12), 3, 4, nrow=4)) 
#           bottom, left, top, right
  par(mar = c(3.2, 3.6, 1.5, 1))   
  for (i in 2:(nmv+1))                 # mv values start in column 2 of mcmcres$mv
   {
# min and max to evaluate for prior
    minmax=c(min(mcmcres$mv[nburnin:niter,i]),max(mcmcres$mv[nburnin:niter,i]))
    fvals=seq(minmax[1],minmax[2],(minmax[2]-minmax[1])/500)
    if(i != iu)
     {
      hist(mcmcres$mv[nburnin:niter,i],freq=FALSE,breaks=100,col="#FF000046",border=NA,xlab="",ylab="",main="")
      box(bty="l")
      mtext(paste0(colnames(mcmcres$mv)[i], " (arcsec/yr)"), side=1,line=2.3,cex=cexset)
      lines(fvals,.comppriorpdf1(mvlabels[i-1],fvals,mvopt,pdfpara['tGa'],pdfpara['umin'],pdfpara['umax']),col="blue",lwd=1.5)
     }else{ 
      hist(100*mcmcres$mv[nburnin:niter,i],freq=FALSE,breaks=100,col="#FF000046",border=NA,xlab="",ylab="",main="")
      box(bty="l")
      mtext(paste0(colnames(mcmcres$mv)[i], " (cm/kyr)"), side=1,line=2.3,cex=cexset)
      lines(100*fvals,.comppriorpdf1(mvlabels[i-1],100*fvals,mvopt,pdfpara['tGa'],100*pdfpara['umin'],100*pdfpara['umax']),col="blue",lwd=1.5)
     }
  }


# PLOT 3a: posterior histograms of astronomical periods
  if(mvopt == 0) 
   {
    gfreqs=mcmcres$mv[nburnin:niter,2:6]
    sfreqs=mcmcres$mv[nburnin:niter,7:11]
    kfreqs=mcmcres$mv[nburnin:niter,ik]
    ap=as.matrix(calcPeriods(g=gfreqs,s=sfreqs,k=kfreqs,opt=2,output=1))
   }
  if(mvopt == 1) 
   {
    gfreqs=mcmcres$mv[nburnin:niter,2:6]
# set s1,s2 and s6 to NA for calcPeriods
    sfreqs0=mcmcres$mv[nburnin:niter,7:8]
    sfreqs=matrix(NA, (niter-nburnin+1), 5)
    sfreqs[,3:4]=sfreqs0[,1:2]
    kfreqs=mcmcres$mv[nburnin:niter,ik]
    ap=as.matrix(calcPeriods(g=gfreqs,s=sfreqs,k=kfreqs,opt=2,output=1))
# remove NA columns output by calcPeriods, associated with s1,s2 and s6 terms
    ap=ap[,-c(6,9,10)]
   }
  if(mvopt == 2)   
   {
    gfreqs=mcmcres$mv[nburnin:niter,2:6]
    kfreqs=mcmcres$mv[nburnin:niter,ik]
    ap=as.matrix(calcPeriods(g=gfreqs,k=kfreqs,opt=2,output=1))
   }

  dev.new(title=c("Posterior PDFs: Milankovitch periods"))
  layoutC <- layout(matrix(c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16), 4, 4, nrow=4)) 
#             bottom, left, top, right
  par(mar = c(3.2, 3.6, 1.5, 1))   
  for (i in 1:dim(ap)[2])   
   {
    hist(ap[,i],freq=FALSE,breaks=100,col="#FF000046",border=NA,xlab="",ylab="",main="")
    box(bty="l")
    mtext(paste0(colnames(ap)[i], " (kyr)"), side=1,line=2.3,cex=cexset)
# NOTE: compare Table 1 Hinnov 2013, Table 1 Malinverno & Meyers 2024, Table 7 (+ Table 5) Laskar et al. 2004
# these values from Malinverno & Meyers 2024
#      if(colnames(ap)[i]=="g2-g5") abline(v=406.2,lty=3,lwd=2,col="blue")
#      if(colnames(ap)[i]=="g3-g2") abline(v=132.0,lty=3,lwd=2,col="blue")
#      if(colnames(ap)[i]=="g4-g2") abline(v=124.0,lty=3,lwd=2,col="blue")
#      if(colnames(ap)[i]=="g3-g5") abline(v=99.6,lty=3,lwd=2,col="blue")
#      if(colnames(ap)[i]=="g4-g5") abline(v=95.0,lty=3,lwd=2,col="blue") 
#      if(colnames(ap)[i]=="k+s6") abline(v=53.7,lty=3,lwd=2,col="blue") 
#      if(colnames(ap)[i]=="k+s3") abline(v=40.9,lty=3,lwd=2,col="blue") 
#      if(colnames(ap)[i]=="k+s4") abline(v=39.6,lty=3,lwd=2,col="blue") 
#      if(colnames(ap)[i]=="k+s2") abline(v=29.6,lty=3,lwd=2,col="blue")     
#      if(colnames(ap)[i]=="k+s1") abline(v=28.9,lty=3,lwd=2,col="blue")     
#      if(colnames(ap)[i]=="k+g5") abline(v=23.7,lty=3,lwd=2,col="blue") 
#      if(colnames(ap)[i]=="k+g1") abline(v=23.0,lty=3,lwd=2,col="blue") 
#      if(colnames(ap)[i]=="k+g2") abline(v=22.4,lty=3,lwd=2,col="blue") 
#      if(colnames(ap)[i]=="k+g3") abline(v=19.1,lty=3,lwd=2,col="blue")  
#      if(colnames(ap)[i]=="k+g4") abline(v=19.0,lty=3,lwd=2,col="blue")  
   }


# PLOT 3b: posterior histograms of grand cycle periods
  if(mvopt == 0 || mvopt == 1) 
   {    
    dev.new(title=c("Posterior PDFs: Grand cycles"))
    layoutD <- layout(matrix(c(1,2), 1, 2, nrow=2)) 
#   hist(g4mg3samplesort,freq=FALSE,breaks=100,col="#FF000046",border=NA,xlab="",ylab="",main="")
# here we are restricting plotting range, to guard against extreme outliers
    hist(g4mg3samplesort[ilow999:ihigh999],freq=FALSE,breaks=100,col="#FF000046",border=NA,xlab="",ylab="",main="",xlim=c(min(0,g4mg3samplesort[ilow999]),max(s4ms3samplesort[ihigh999],g4mg3samplesort[ihigh999])) )
    box(bty="l") 
    mtext("g4-g3 (Myr)", side=1,line=2.3,cex=cexset)
#   hist(s4ms3samplesort,freq=FALSE,breaks=100,col="#FF000046",border=NA,xlab="",ylab="",main="")
# here we are restricting plotting range, to guard against extreme outliers
    hist(s4ms3samplesort[ilow999:ihigh999],freq=FALSE,breaks=100,col="#FF000046",border=NA,xlab="",ylab="",main="",xlim=c(min(0,s4ms3samplesort[ilow999]),max(s4ms3samplesort[ihigh999],g4mg3samplesort[ihigh999])) )
    box(bty="l") 
    mtext("s4-s3 (Myr)", side=1,line=2.3,cex=cexset)   
   } 

  if(mvopt == 2) 
   {    
    dev.new(title=c("Posterior PDF: Grand cycle"))
#   hist(g4mg3samplesort,freq=FALSE,breaks=100,col="#FF000046",border=NA,xlab="",ylab="",main="")
# here we are restricting plotting range, to guard against extreme outliers
    hist(g4mg3samplesort[ilow999:ihigh999],freq=FALSE,breaks=100,col="#FF000046",border=NA,xlab="",ylab="",main="",xlim=c(min(0,g4mg3samplesort[ilow999]),g4mg3samplesort[ihigh999]) )
    box(bty="l") 
    mtext("g4-g3 (Myr)", side=1,line=2.3,cex=cexset)
   } 
    
    
# PLOT 4: crossplots
  dev.new(title=c("Correlations"))
  panel.cor <- function(x, y, digits = 2, ...)
   {
    par(usr = c(0, 1, 0, 1))
    r <- cor(x, y)
    if(r>0) color="red"
    if(r<0) color="blue"
    txt <- format(c(r, 0.123456789), digits = digits)[1] 
    text(0.5, 0.5, txt,cex=0.9,font=4,col=color)
   }
  pairs(mcmcres$mv[seq(nburnin,niter,by=ceiling((niter-nburnin)/2000)),2:(nmv+1)],upper.panel=panel.cor,gap=0.5,col="#0000001E",pch=19,cex=.3,cex.labels=1.2,xaxt="n",yaxt="n",row1attop=T)
  
  
  
# PLOT 5: plot fit to the data for MAP sampled value of mv
# note that in mcmcres$mv, the mvs start with second column  
  .plotmvdatafit(dat,fc_dat,mvmap[2:(nmv+1)],pdfpara["mvopt"],iu-1,pdfpara['roll'],pdfpara['deltafp'],pdfpara['nolikfdenv'])

# end function timeOptBMCMCplot
 }