### This function is a component of astrochron: An R Package for Astrochronology
### Copyright (C) 2018 Stephen R. Meyers
###
###########################################################################
### multiTest function - (SRM: November 19-30, 2017; December 4, 2017;
###                            April 13, 2018; April 17, 2018; July 7, 2018;
###                            January 8, 2019)
###
### apply multiple comparisons correction to spectrum confidence levels
###########################################################################

multiTest <- function (spec,flow=NULL,fhigh=NULL,pl=T,output=T,genplot=T,verbose=T)
{

# spec should have two columns: frequency and confidence level
if(verbose) cat("\n----- Adjusting spectral p-values for multiple comparisons -----\n")

spec <- data.frame(spec)

# if there are 2 columns, they must be in the order: frequency, uncorrected confidence level
if(length(spec)==2)
 {
  freq <- spec[,1]
  rawPvals = (100-spec[,2])/100
 }

# if there are 8 columns, we assume the results come from lowspec,mtm,mtmPL,or mtmML96
# note that periodgram has a different order!
if(length(spec)==8)
 {
   freq <- spec[,1]
   rawPvals = (100-spec[,4])/100
 }
 
# if there are 9 columns, we assume the results come from periodogram
if(length(spec)==9)
 {
   freq <- spec[,1]
   rawPvals = (100-spec[,5])/100
 }

if(length(spec)<2 || length(spec)>9)
 {
    cat("\n**** ERROR: spec does not have the correct column structure.\n")
    stop("**** TERMINATING NOW!")
 }

if(length(spec)>2 && length(spec)<8)
 {
    cat("\n**** ERROR: spec does not have the correct column structure.\n")
    stop("**** TERMINATING NOW!")
 }

ifreq=length(freq)

if(is.null(flow)) 
 {
   flow=freq[1]
   if(verbose) cat("* Using default lower frequency bound=", flow,"\n")
 }

if(is.null(fhigh)) 
 {
   fhigh=freq[ifreq]
   if(verbose) cat("* Using default upper frequency bound=", fhigh,"\n")
 }
ifreq2=0
for (i in 1:length(flow)) ifreq2=append(ifreq2,which( (freq >= flow[i]) & (freq <= fhigh[i]) ))
# remove first entry, which is 0
ifreq2=ifreq2[2:length(ifreq2)]
# remove duplicate values
ifreq2=unique(ifreq2)
pvals=rawPvals[ifreq2]

# fdr is the same as bh
fdr=p.adjust(pvals,method=c("fdr"))
BY=p.adjust(pvals,method=c("BY"))
hommel=p.adjust(pvals,method=c("hommel"))
hochberg=p.adjust(pvals,method=c("hochberg"))
holm=p.adjust(pvals,method=c("holm"))
bonferroni=p.adjust(pvals,method=c("bonferroni"))


out = data.frame(cbind(freq[ifreq2],1/freq[ifreq2],rawPvals[ifreq2],fdr,BY,hommel,hochberg,holm,bonferroni))
colnames(out) = c("Frequency","Period","Uncorrected p-value","FDR-BH","FDR-BY","Hommel","Hochberg","Holm","Bonferroni")

plSwitch=T
if(min(flow)==freq[1] && max(fhigh)==freq[ifreq]) plSwitch=F

setLwd=0
if(pl) setLwd=1.5

if(genplot)
 {
   par(mfrow=c(3,2))
   plot(freq,rawPvals,type="l",lwd=setLwd,col="gray",xlab="Frequency",ylab="FDR-BH",main="(a) False Discovery Rate-BH",bty="n",ylim=c(0,1))
#   if(plSwitch) rect(flow,-1,fhigh,max(out[,5]),col="#BEBEBE50",border="NA")
   points(out[,1],out[,4],pch=20,col="red")
   abline(h=0.05,lty=3,col="blue")
   text(max(freq),0.05,"0.05",font=2,pos=1,offset=0.15,col="blue")
   abline(h=0.1,lty=3)
   text(max(freq),0.1,"0.1",font=2,pos=3,offset=0.1)

   plot(freq,rawPvals,type="l",lwd=setLwd,col="gray",xlab="Frequency",ylab="FDR-BY",main="(b) False Discovery Rate-BY",bty="n",ylim=c(0,1))
#   if(plSwitch) rect(flow,-1,fhigh,max(out[,5]),col="#BEBEBE50",border="NA")
   points(out[,1],out[,5],pch=20,col="red")
   abline(h=0.05,lty=3,col="blue")
   text(max(freq),0.05,"0.05",font=2,pos=1,offset=0.15,col="blue")
   abline(h=0.1,lty=3)
   text(max(freq),0.1,"0.1",font=2,pos=3,offset=0.1)

   plot(freq,rawPvals,type="l",lwd=setLwd,col="gray",xlab="Frequency",ylab="Hommel-corrected",main="(c) Hommel-corrected p-values",bty="n",ylim=c(0,1))
#   if(plSwitch) rect(flow,-1,fhigh,max(rawPvals),col="#BEBEBE50",border="NA")
   points(out[,1],out[,6],pch=20,col="red")
   abline(h=0.05,lty=3,col="blue")
   text(max(freq),0.05,"0.05",font=2,pos=1,offset=0.15,col="blue")
   abline(h=0.1,lty=3)
   text(max(freq),0.1,"0.1",font=2,pos=3,offset=0.1)

   plot(freq,rawPvals,type="l",lwd=setLwd,col="gray",xlab="Frequency",ylab="Hochberg-corrected",main="(d) Hochberg-corrected p-values",bty="n",ylim=c(0,1))
#   if(plSwitch) rect(flow,-1,fhigh,max(out[,6]),col="#BEBEBE50",border="NA")
   points(out[,1],out[,7],pch=20,col="red")
   abline(h=0.05,lty=3,col="blue")
   text(max(freq),0.05,"0.05",font=2,pos=1,offset=0.15,col="blue")
   abline(h=0.1,lty=3)
   text(max(freq),0.1,"0.1",font=2,pos=3,offset=0.1)

   plot(freq,rawPvals,type="l",lwd=setLwd,col="gray",xlab="Frequency",ylab="Holm-corrected",main="(e) Holm-corrected p-values",bty="n",ylim=c(0,1))
#   if(plSwitch) rect(flow,-1,fhigh,max(out[,7]),col="#BEBEBE50",border="NA")
   points(out[,1],out[,8],pch=20,col="red")
   abline(h=0.05,lty=3,col="blue")
   text(max(freq),0.05,"0.05",font=2,pos=1,offset=0.15,col="blue")
   abline(h=0.1,lty=3)
   text(max(freq),0.1,"0.1",font=2,pos=3,offset=0.1)

   plot(freq,rawPvals,type="l",lwd=setLwd,col="gray",xlab="Frequency",ylab="Bonferroni-corrected",main="(f) Bonferroni-corrected p-values",bty="n",ylim=c(0,1))
#   if(plSwitch) rect(flow,-1,fhigh,max(out[,8]),col="#BEBEBE50",border="NA")
   points(out[,1],out[,9],pch=20,col="red")
   abline(h=0.05,lty=3,col="blue")
   text(max(freq),0.05,"0.05",font=2,pos=1,offset=0.15,col="blue")
   abline(h=0.1,lty=3)
   text(max(freq),0.1,"0.1",font=2,pos=3,offset=0.1)
 }
  
if(output) return(out)
#### END function multiTest
}
