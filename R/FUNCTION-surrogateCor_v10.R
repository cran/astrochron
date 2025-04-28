### This function is a component of astrochron: An R Package for Astrochronology
### Copyright (C) 2023 Stephen R. Meyers
###
###########################################################################
### function surrogateCor: calculate signifiance of correlation coefficient
###       for serially correlated data using method of Ebisuzaki W. (1997);
###       A method to estimate the statistical significance of a correlation 
###       when the data are serially correlated. J Climate, 10, 2147-2153).
###       This version allows the comparision of variables with
###       two different sampling grids.
###                            (SRM: April 5-8, 2014; May 1, 2014; May 2, 2014;
###                                   October 10, 2014; November 13, 2015; 
###                                   Dec. 25, 2015; May 20-25, 2016; 
###                                   June 7, 2016; July 1, 2016; July 22, 2016;
###                                   January 9, 2019; January 14, 2021;
###                                   April 12, 2023; November 10, 2024)
###
### see also http://www.mathworks.com/matlabcentral/fileexchange/10881-weaclim/content/ebisuzaki.m
###########################################################################



surrogateCor <- function (dat1=NULL,dat2=NULL,firstDiff=F,cormethod=1,nsim=1000,output=2,genplot=T,check=T,verbose=T)
{

if(verbose) 
 {
   cat("\n--------- CALCULATING CORRELATION COEFFICIENT ---------\n")
   cat("     AND SIGNIFICANCE FOR SERIALLY CORRELATED DATA\n")
 }  

if(is.null(dat1) || is.null(dat2)) { stop("\n**** ERROR: dat1 and dat2 must be provided. TERMINATING NOW!") }

dat1<- data.frame(dat1)
dat2<- data.frame(dat2)

if(check)
 {
# check if dat1 and dat2 have same number of columns
   if(length(dat1)!= length(dat2))
    {
      cat("\n**** ERROR: dat1 and dat2 have different number of columns\n")
      stop("**** TERMINATING NOW!")
    }
 }  

############################################
# if both dat1 and dat2 have two columns
############################################
if(length(dat1)==2 && length(dat2)==2)
 {
if(verbose) cat("\n * IMPLEMENTING RESAMPLING ALGORITHM\n")
# sort to ensure increasing depth/height/time
if(check)
 {
   if(verbose) cat("\n * Sorting datasets to ensure increasing order, removing empty entries\n")
   dat1 <- dat1[order(dat1[,1],na.last=NA,decreasing=F),]
   dat2 <- dat2[order(dat2[,1],na.last=NA,decreasing=F),]
 } 

n1=length(dat1[,1])
n2=length(dat2[,1])
if(verbose) cat(" * Number of data points in dat1=",n1,"\n")
if(verbose) cat(" * Number of data points in dat2=",n2,"\n")

if(check)
 {
# check for duplicate depths/heights in dat1
   dx1=dat1[2:n1,1]-dat1[1:(n1-1),1]
   if(min(dx1) == 0)
     {
       cat("\n**** ERROR: duplicate depth/height datum found in dat1\n")
       stop("**** TERMINATING NOW!")
     }  

### check for duplicate depths/heights in dat2
   dx2=dat2[2:n2,1]-dat2[1:(n2-1),1]
   if(min(dx2) == 0)
     {
       cat("\n**** ERROR: duplicate depth/height datum found in dat2\n")
       stop("**** TERMINATING NOW!")
     }  
# end check
 }

# isolate data from interval that is shared by dat1 and dat2
xmin= max( min(dat1[,1]), min(dat2[,1]) )
xmax= min( max(dat1[,1]), max(dat2[,1]) )
if(verbose) cat(" * Isolating data between the shared interval of", xmin, "and", xmax, "\n")
dat1S= subset( dat1, (dat1[1] >= xmin) & (dat1[1] <= xmax) )
dat2S= subset( dat2, (dat2[1] >= xmin) & (dat2[1] <= xmax) )
n1S= length(dat1S[,1])
n2S= length(dat2S[,1])
if(verbose) cat(" * New number of data points in dat1=",n1S,"\n")
if(verbose) cat(" * New number of data points in dat2=",n2S,"\n")

# ensure that values just beyond common interval are available for interpolation purposes.
n1Snew= n1S
n2Snew= n2S

if(n1S <= n2S) 
 {
   if(dat1S[1,1] < dat2S[1,1])
     {
       ii=which(dat2[1]==dat2S[1,1])
       a=append(dat2S[,1],dat2[ii-1,1], after=0)
       b=append(dat2S[,2],dat2[ii-1,2], after=0)
       dat2S=data.frame(cbind(a,b))
       n2Snew= length(dat2S[,1])
     }
     
   if(dat1S[n1S,1] > dat2S[n2S,1])
     {
       ii=which(dat2[1]==dat2S[n2Snew,1])
       a=append(dat2S[,1],dat2[ii+1,1], after=n2Snew)
       b=append(dat2S[,2],dat2[ii+1,2], after=n2Snew)
       dat2S=data.frame(cbind(a,b))
       n2Snew= length(dat2S[,1])
     }     
 }
if(n1S > n2S) 
 {
   if(dat1S[n1S,1] < dat2S[n2S,1])
     {
       ii=which(dat1[1]==dat1S[n1Snew,1])
       a=append(dat1S[,1],dat1[ii+1,1], after=n1Snew)
       b=append(dat1S[,2],dat1[ii+1,2], after=n1Snew)
       dat1S=data.frame(cbind(a,b))
       n1Snew= length(dat1S[,1])
     } 
   if(dat1S[1,1] > dat2S[1,1])
     {
       ii=which(dat1[1]==dat1S[1,1])
       a=append(dat1S[,1],dat1[ii-1,1], after=0)
       b=append(dat1S[,2],dat1[ii-1,2], after=0)
       dat1S=data.frame(cbind(a,b))
       n1Snew= length(dat2S[,1])
     }
 }   
              
x1S=dat1S[,2]
x2S=dat2S[,2]
# remove mean
x1S=x1S-mean(x1S)
x2S=x2S-mean(x2S)

# resample on grid from series with fewer data points
if(n1S <= n2S) 
 {
   if(verbose) cat("\n * Resampling dat2 on dat1 sample grid.\n")
   dat2SR=resample(dat2S,dat1S[,1],verbose=F,check=F,genplot=F)
   dat1SR=dat1S
   denTitle1=c("Distribution of Values for dat1")
   denTitle2=c("Distribution of Resampled Values for dat2")
   stratTitle1=c("Stratigraphic Series for dat1")
   stratTitle2=c("Stratigraphic Series for dat2 (red=resampled)")
   plSwitch=2
 } 
if(n1S > n2S) 
 {
   if(verbose) cat(" * Resampling dat1 on dat2 sample grid.\n")
   dat1SR=resample(dat1S,dat2S[,1],verbose=F,check=F,genplot=F)
   dat2SR=dat2S
   denTitle1=c("Distribution of Resampled Values for dat1")
   denTitle2=c("Distribution of Values for dat2")
   plSwitch=1
 }  

x1SR=dat1SR[,2]
x2SR=dat2SR[,2]
}

############################################
# if dat1 and dat2 have one column
############################################
if(length(dat1)==1 && length(dat2)==1)
{
 if(verbose) cat("\n * NO RESAMPLING APPLIED\n")
 n1=length(dat1[,1])
 n2=length(dat2[,1])
 
if(check)
  {
 if(n1 != n2)
      {
       cat("\n**** ERROR: dat1 and dat2 must have the same number of points")
       stop("**** TERMINATING NOW!")
     }  
  }
  
 if(verbose) cat("\n * Number of data points per variable=",n1,"\n")

 x1SR=dat1[,1]
 x2SR=dat2[,1]
 plSwitch=3
 denTitle1=c("Distribution of Values for dat1")
 denTitle2=c("Distribution of Values for dat2")
}

# remove mean
x1SR=x1SR-mean(x1SR)
x2SR=x2SR-mean(x2SR)

if(length(dat1)==1 && length(dat2)==1)
 {
  x1S=x1SR
  x2S=x2SR
 }
  
# caculate first differences if desired
if(firstDiff) 
 {
  if(verbose) cat("\n * Calculating first differences.\n")
  x1SR=diff(x1SR)
  x2SR=diff(x2SR)
  }

if(cormethod == 1) method = c("pearson")
if(cormethod == 2) method = c("spearman")
if(cormethod == 3) method = c("kendall")

# calculate correlation coefficient for data
datcor=cor(x1SR,x2SR,method=method)
if(verbose) cat("\n * Data correlation coefficient =",datcor,"\n")


if(nsim >1)
{
### Simulation
# generate surrogates
if(verbose) cat("\n * Generating phase randomized surrogates for variable 1\n")
x1sur <- surrogates(x1S,nsim=nsim,verbose=F)
if(verbose) cat(" * Generating phase randomized surrogates for variable 2\n")
x2sur <- surrogates(x2S,nsim=nsim,verbose=F)

if(verbose) cat(" * Calculating correlation coefficient for surrogates\n")
surrcor<-double(nsim)
# calculate correlation coefficient for simulations
for (i in 1:nsim)
 {

  if(length(dat1)==2 && length(dat2)==2)
   {
# resample on grid from series with fewer data points
    if(n1S <= n2S) 
      {
        x1surNow=x1sur[,i]
        x2surNow=resample(cbind(dat2S[,1],x2sur[,i]),dat1S[,1],verbose=F,check=F,genplot=F)[2]
        x2surNow=x2surNow[,1] 
      } 
    if(n1S > n2S) 
      {
        x1surNow=resample(cbind(dat1S[,1],x1sur[,i]),dat2S[,1],verbose=F,check=F,genplot=F)[2]
        x1surNow=x1surNow[,1]
        x2surNow=x2sur[,i]
      }
   }

  if(length(dat1)==1 && length(dat2)==1)
   {
     x1surNow=x1sur[,i]
     x2surNow=x2sur[,i]
   }
  
    if(firstDiff) 
      {
        x1surNow=diff(x1surNow)
        x2surNow=diff(x2surNow)
      }
     
    surrcor[i]=cor(x1surNow,x2surNow,method=method)
 }  

if(verbose) cat(" * Evaluating significance level\n")
#   number of results with abs(surrcor) > abs(datcor)
    pnumgt = sum(abs(surrcor)>abs(datcor))    
if(verbose) cat("\n * Number of simulations with |r<sim>| > |r<dat>| = ", pnumgt,"\n")
    pvalue=pnumgt/nsim
    if(pvalue >= (10/nsim) && (10/nsim) <=1 ) if(verbose) cat(" * P-value =", pvalue,"\n")  
    if((10/nsim) > 1 ) 
      {
       pvalue=1
       if(verbose) cat(" * P-value = 1 \n")  
      } 
    if(pvalue < (10/nsim) && (10/nsim) <=1 ) 
     {
       pvalue=10/nsim
       if(verbose) cat(" * P-value <", pvalue,"\n")
     }  
     
# end if nsim >1
}

if(genplot)
 {
    par(mfrow=c(3,2))
    if(plSwitch==1)
     {
       plot(dat1S,col="black",type="b",xlim=c(xmin,xmax),xlab="Location",ylab="dat1 data value",main="Stratigraphic Series for dat1 (red=resampled)"); points(dat1SR,col="red",cex=0.7)
       plot(dat2S,col="black",type="b",xlim=c(xmin,xmax),xlab="Location",ylab="dat2 data value",main="Stratigraphic Series for dat2")
     }
    if(plSwitch==2)
     {
       plot(dat1S,col="black",type="b",xlim=c(xmin,xmax),xlab="Location",ylab="dat1 data value",main="Stratigraphic Series for dat1")
       plot(dat2S,col="black",type="b",xlim=c(xmin,xmax),xlab="Location",ylab="dat2 data value",main="Stratigraphic Series for dat2 (red=resampled)"); points(dat2SR,col="red",cex=0.7)
     }

    if(plSwitch==3)
     {
       plot(dat1[,1],col="black",type="b",xlab="Location",ylab="dat1 data value",main="Stratigraphic Series for dat1")
       plot(dat2[,1],col="black",type="b",xlab="Location",ylab="dat2 data value",main="Stratigraphic Series for dat2")
     }

    if(!firstDiff)
     {
       if(plSwitch !=3)
        { 
          hist(dat1SR[,2],freq=F, xlab="dat1 data value",main=denTitle1); lines(density(dat1SR[,2], bw="nrd0"),col="red",lwd=1.5)
          hist(dat2SR[,2],freq=F, xlab="dat2 data value",main=denTitle2); lines(density(dat2SR[,2], bw="nrd0"),col="red",lwd=1.5)
          plot(dat1SR[,2],dat2SR[,2], xlab="dat1 data value", ylab="dat2 data value",main="dat1 vs. dat2")  
        }         
       if(plSwitch ==3)
        { 
          hist(dat1[,1],freq=F, xlab="dat1 data value",main=denTitle1); lines(density(dat1[,1], bw="nrd0"),col="red",lwd=1.5)
          hist(dat2[,1],freq=F, xlab="dat2 data value",main=denTitle2); lines(density(dat2[,1], bw="nrd0"),col="red",lwd=1.5)
          plot(dat1[,1],dat2[,1], xlab="dat1 data value", ylab="dat2 data value",main="dat1 vs. dat2")  
        }         
     } 
    if(firstDiff) 
     {
       hist(x1SR,freq=F, xlab="dat1 data value",main=denTitle1); lines(density(x1SR, bw="nrd0"),col="red",lwd=1.5)
       hist(x2SR,freq=F, xlab="dat2 data value",main=denTitle2); lines(density(x2SR, bw="nrd0"),col="red",lwd=1.5)
       plot(x1SR,x2SR, xlab="dat1 data value", ylab="dat2 data value",main="dat1 vs. dat2 (first differences)")
     }

    if(nsim>1)
     {
       plotMin=min(datcor,surrcor)
       plotMax=max(datcor,surrcor)
       denRes=density(surrcor, bw="nrd0")
       if(datcor < 0) plotPos=4
       if(datcor >= 0) plotPos=2
# set up plot       
       plot(denRes,col="red",xlim=c(plotMin,plotMax),xlab="Simulated Correlation Coefficient",main="Monte Carlo Simulation Results",lwd=1.5)
# fill
       polygon(denRes,col="gray",border=NA)
# replot line to show full widith
       lines(denRes,col="red",lwd=1.5)
       abline(v=datcor,lwd=2,lty=6,col="blue")
       text(datcor,0.6*max(denRes$y),labels=round(datcor,digits=3),pos=plotPos,cex=1.4,col="blue")
     }
  }
  
  if(output==1 && nsim>1) return(surrcor)
  if(output==2)
   {
     if(nsim>1) out= data.frame( cbind(datcor,pvalue) )
     if(nsim<=1) out= data.frame(datcor)
     return(out)
   }
  if(output==3 && length(dat1)==2 && length(dat2)==2) 
   {
     if(!firstDiff)
      {
        out= data.frame( cbind(dat1SR[,1],dat1SR[,2],dat2SR[,2]) )
        colnames(out) <- c("location","dat1","dat2")
      }
    if(firstDiff)
      {
        out= data.frame( cbind(dat1SR[2:length(dat1SR[,1]),1],x1SR,x2SR) )
        colnames(out) <- c("location","dat1_diff","dat2_diff")
      }     
     return(out)
   }
   
# end function surrogateCor
}