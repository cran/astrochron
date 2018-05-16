### This function is a component of astrochron: An R Package for Astrochronology
### Copyright (C) 2017 Stephen R. Meyers
###
###########################################################################
### slideCor function - (SRM: June 28, 2013; May 5-8, 2017; May 20-26, 2017;
###                           June 12, 2017; August 31, 2017; October 11, 2017)
###
### given two data series, slide one data series past the other
### to identify the shift required for maximum correlation
###########################################################################

slideCor <- function (dat1,dat2,rev=F,cormethod=1,minpts=5,output=T,genplot=T,verbose=T)
{

if(verbose) 
 {
   cat("\n----- IDENTIFY OPTIMAL TEMPORAL OR SPATIAL SHIFT NECESSARY -----\n")
   cat("          TO MAXIMIZE CORRELATION BETWEEN TWO DATA SERIES\n")
 }  

dat1 <- data.frame(dat1)
npts1 <- length(dat1[,1]) 
dt1 <- dat1[2,1]-dat1[1,1]

if(dt1<1)
 { 
    if (verbose) cat("\n * Sorting stratigraphic series 1 into increasing height/depth/time, removing empty entries\n")
    dat1 <- dat1[order(dat1[1], na.last = NA, decreasing = F), ]
    dt1 <- dat1[2,1]-dat1[1,1]
    npts1 <- as.integer( length(dat1[,1]) )
 }

# error checking 
dtest <- dat1[2:npts1,1]-dat1[1:(npts1-1),1] 
epsm=1e-9
if( (max(dtest)-min(dtest)) > epsm ) 
 {
   cat("\n**** ERROR: dat1 sampling interval is not uniform.\n")
   stop("**** TERMINATING NOW!")
 }

if (verbose) 
 {
   cat("\n * Number of data points in stratigraphic series 1:",npts1,"\n")
   cat(" * Stratigraphic series length (space or time):",(npts1-1)*dt1,"\n")
   cat(" * Sampling interval (space or time):",dt1,"\n")
 }

dat2 <- data.frame(dat2)
npts2 <- length(dat2[,1]) 
dt2 <- dat2[2,1]-dat2[1,1]

if(rev) dat2[2]= -1 * dat2[2]

if(dt2<1)
 { 
    if (verbose) cat("\n * Sorting stratigraphic series 2 into increasing height/depth/time, removing empty entries\n")
    dat2 <- dat2[order(dat2[1], na.last = NA, decreasing = F), ]
    dt2 <- dat2[2,1]-dat2[1,1]
    npts2 <- as.integer( length(dat2[,1]) )
 }
 
dtest <- dat2[2:npts2,1]-dat2[1:(npts2-1),1] 
epsm=1e-9
if( (max(dtest)-min(dtest)) > epsm ) 
 {
       cat("\n**** ERROR: dat2 sampling interval is not uniform.\n")
       stop("**** TERMINATING NOW!")
 }

if (verbose) 
 {
   cat("\n * Number of data points in stratigraphic series 2:",npts2,"\n")
   cat(" * Stratigraphic series length (space or time):",(npts2-1)*dt2,"\n")
   cat(" * Sampling interval (space or time):",dt2,"\n")
 }

if(dt1 == dt2) dtOut=dt1 
# resample to coarsest resolution
if(dt1 != dt2) 
  { 
    if(dt1 > dt2) 
      {
        if(verbose) cat("\n * Decimating series 2 to sampling interval of series 1\n")
        dat2 <- linterp(dat2,dt=dt1,genplot=F,verbose=F)
        npts2 <- length(dat2[,1]) 
        dtOut <- dt1
       } 
    if(dt1 < dt2) 
      {
        if(verbose) cat(" * Decimating series 1 to sampling interval of series 2\n")
        dat1 <- linterp(dat1,dt=dt2,genplot=F,verbose=F)  
        npts1 <- length(dat1[,1]) 
        dtOut <- dt2
       } 
  }

if(npts1 >= npts2)
 {
   datL=dat1
   datS=dat2
   nptsL=npts1
   nptsS=npts2
 }  

if(npts1 < npts2)
 {
   datS=dat1
   datL=dat2
   nptsS=npts1
   nptsL=npts2
 }  

dpts=nptsL-nptsS

# more error checking
if(minpts<2)
 {
       cat("\n**** WARNING: minpts must be set to a value >= 2. Resetting to 2.\n")
       minpts=2
 }

if(minpts>nptsS)
 {
       cat("\n**** WARNING: minpts must be set to a value <=",nptsS, "  Resetting to 2.\n")
       minpts=2
 }

if(cormethod == 1) method=c("pearson")
if(cormethod == 2) method=c("spearman")
if(cormethod == 3) method=c("kendall")

ptsOut=nptsL+nptsS-3
res <- double(ptsOut)
shift <- double(ptsOut)
pts <- double(ptsOut)

# result ticker
k = 1

kTrunc = 1
for (i in 1:ptsOut)
 {
# ticker for short series
    j=nptsS-k
# ticker for long series
    jj=k+1        
    if(j > 0) 
        {
           a=datS[j:nptsS,2]
           b=datL[1:jj,2]
           pts[k] = nptsS-j+1
        }

   if(j <= 0 && jj <= nptsL) 
        {
           a=datS[1:nptsS,2]
           b=datL[((-1*j)+2):jj,2]
           pts[k] = nptsS  
        }
        
   if(j <= 0 && jj > nptsL) 
       {
           a=datS[1:(nptsS-kTrunc),2]  
           b=datL[(nptsL-nptsS+kTrunc+1):nptsL,2]
           pts[k] = nptsS-kTrunc 
           kTrunc = kTrunc + 1
       }
         
   res[k] = cor(a,b,method=method)
   shift[k] = (k - nptsS + 1)*dtOut
   k = k + 1
 }

out = data.frame(cbind(shift,res,pts))
out = subset(out, (out[3] >= minpts )) 
colnames(out) <- c("shift","correlation_coefficient","points")

# identify maximum r-value
i1=which.max(out[,2])
# identify maximum r2-value
i2=which.max(out[,2]^2)

# check for multiple maxima
resMax=which(out[,2]==out[i1,2])
resMax2=which(out[,2]^2==out[i2,2]^2)

sumMax = sum( out[,2]==out[i1,2])
if(length(resMax)>1 && verbose) cat("\n**** WARNING:", sumMax,"maxima detected in r-values. The first occurence will be used.\n") 
sumMax2 = sum( out[,2]^2==out[i2,2]^2 )
if(length(resMax2)>1 && verbose) cat("\n**** WARNING:", sumMax2,"maxima detected in r2-values.\n") 

if(verbose)
 {
   cat("\n * Optimal shift identified using r-values ( r=",out[i1,2],"): \n")
   print(out[resMax,1])
   cat(" * Optimal shift identified using r2-values ( r2=",out[i2,2]^2,"): \n")
   print(out[resMax2,1])
 }
 
if(genplot)
 {
     dev.new(height = 7.7, width = 5.7, title = "slideCor Results")
     par(mfrow=c(5,1),mar=c(3.1, 4.1, 4.1, 2.1))
     plot(dat1, cex=.5,xlab="",ylab="",main="Data Series 1",bty="n")
     lines(dat1)
     mtext("Location (spatial or temporal units)",side=1,line=2,cex=0.7)
     mtext("Value",side=2,line=2,cex=0.7)
     par(mar=c(4.1, 4.1, 4.1, 2.1))
     plot(dat2, cex=.5,xlab="",ylab="",main="Data Series 2",bty="n")
     lines(dat2,col="red")
     mtext("Location (spatial or temporal units)",side=1,line=2,cex=0.7)
     mtext("Value",side=2,line=2,cex=0.7)
     par(mar=c(4.1, 4.1, 3.1, 2.1))
     plot(datL, cex=.5,xlab="",ylab="",main="Optimal Match",bty="n",xlim=c(min(datL[,1]),max(datL[,1])))

     start=datL[1,1]+out[resMax[1],1]
     end=(datL[1,1]+out[resMax[1],1])+(nptsS-1)*dtOut
     par(new=T)
     plot(seq(from=start,to=end,by=dtOut),datS[,2],type="l",col="red",xlab="",ylab="",main="",bty="n",xaxt="n",yaxt="n",xlim=c(min(datL[,1]),max(datL[,1])))
     mtext("Location (spatial or temporal units)",side=1,line=2,cex=0.7)
     mtext("Value",side=2,line=2,cex=0.7)
     plot(out[,1], out[,2], cex=.5,xlab="",ylab="",main="Sliding Correlation Result",bty="n")
     lines(out[,1], out[,2])
     abline(v=out[resMax,1],col="red",lty=2)
     mtext("Shift (spatial or temporal units)",side=1,line=2,cex=0.7)
     mtext("Correlation",side=2,line=3,cex=0.7)
     mtext("Coefficient",side=2,line=2,cex=0.7)
     plot(out[,1], out[,3], cex=.5,xlab="",ylab="",main="# Data Points Compared",bty="n")  
     mtext("Shift (spatial or temporal units)",side=1,line=2,cex=0.7)
     mtext("# Data Points",side=2,line=2,cex=0.7)
 }   
  
if(output) return(out)

### END function slideCor
}
