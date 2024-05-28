### This function is a component of astrochron: An R Package for Astrochronology
### Copyright (C) 2024 Stephen R. Meyers
###
###########################################################################
### function insoMap : calculate insolation map using Laskar et al. (2004)
###                    solution, with palinsol (SRM: November 7-9, 2022;
###                    July 22-24, 2023; March 12, 2024)
###
###########################################################################


insoMap <- function (t=0,pick=0,S0=1365,output=F,verbose=T) 
{
   if(verbose) cat("\n----- CALCULATE INSOLATION MAP USING LASKAR ET AL. (2004) SOLUTION -----\n")

# error checking
# restrict to 21 million years into future (>= -21000). 
   if(t < (-21000))
    {
     if(verbose) cat("\n**** WARNING: This solution can only be calculated to 21000 kiloyears in the future. Resetting t to -21000.\n")
     t=-21000
    }   

# restrict to 51 million years into past (<= 51000).
   if(t > 51000)
    {
     if(verbose) cat("\n**** WARNING: This solution can only be calculated to 51000 kiloyears in the past. Resetting t to 51000.\n")
     t=51000
    }    

   la2=la04(t*1000)
   M <- Milankovitch(la2)
   obl=la2[1]*180/pi
   ecc=la2[2]
   varpi=la2[3]
   plot(M, plot=contour,main=paste(t,"ka"))
   if(pick!=2) mtext(paste("Tilt=",round(obl,2),"    Eccentricity=",round(ecc,2),"    Perihelion=",round(varpi,2)),side=3,line=0.5,font=2)

# below adapted from palinsol's Insol.R
   lat=attr(M,"lat")
   long=attr(M,"long")
   Col  = c(which(long >= (360-80)) , which(long < (360-80)))
   Col = c(Col, Col[1])
   long2 = c(long, long[1]+360)
   MM = M[ Col,]
# above adapted from palinsol's Insol.R


## this script modified from '?identify' in R
identifyPch <- function(x, y=NULL, n=length(x), pch=19, ...)
{
    xy <- xy.coords(x, y); x <- xy$x; y <- xy$y
    sel <- rep(FALSE, length(x)); res <- integer(0)
    while(sum(sel) < n) {
        ans <- identify(x[!sel], y[!sel], n=1, plot=F, ...)
        if(!length(ans)) break
        ans <- which(!sel)[ans]
        points(x[ans], y[ans], pch = pch, col="black")
        sel[ans] <- TRUE
        res <- c(res, ans)
    }
    res
}

   if(pick==1)
    {
      par(new=T)
      plot(rep(370,length(lat)),lat,xlim=c(0,357.5),ylim=c(-90,90),xaxs="i",yaxs="i",yaxt='n',xaxt='n',bty='n',ylab="",xlab="",pch=1,cex=.5,col="black",xpd = TRUE)
      cat("\n *****  Select by clicking on a point at the right of the graph  *****\n")
      pts <- identifyPch(rep(370,length(lat)),lat,n=1)
      abline(h=lat[pts],col="red",lwd=2,lty=4)
      if(verbose) cat("\n * Extracting record from=",lat[pts],"degrees\n")     
      dev.new()
      plot(long2, MM[,pts],type="l",lwd=2,col="red",xaxt='n',xlab="Month",ylab="",main=paste("Latitude=",round(lat[pts],2),"degrees"))
      mtext(expression(paste("Incoming Solar Radiation W/m"^"2")),side=2,line=2.5)
      mtext(paste("Tilt=",round(obl,2),"    Eccentricity=",round(ecc,2),"    Perihelion=",round(varpi,2)),side=3,font=2)
      axis (1, at=seq(0,12)*30, labels=rep('',13))
      axis (1, at=seq(0,11)*30+15, labels=c('J','F','M','A','M','J','J','A','S','O','N','D'), tick=FALSE)
     }

   if(pick==2)
    {
      par(new=T)
      plot(long2,rep(96,length(long2)),xlim=c(0,357.5),ylim=c(-90,90),xaxs="i",yaxs="i",yaxt='n',xaxt='n',bty='n',ylab="",xlab="",pch=1,cex=.5,col="black",xpd = TRUE)
      cat("\n *****  Select by clicking on a point at the top of the graph  *****\n")
      pts <- identifyPch(long2,rep(96,length(long2)),n=1)
      abline(v=long2[pts],col="red",lwd=2,lty=4)
      if(verbose) cat("\n * Extracting record from=",long2[pts],"degrees\n")     
      dev.new()
      plot(MM[pts,],lat,type="l",lwd=2,col="red",ylab="Latitude (degrees)",xlab=expression(paste("Incoming Solar Radiation W/m"^"2")),main=paste("True Solar longitude=",round(long2[pts],2)))
      mtext(paste("Tilt=",round(obl,2),"    Eccentricity=",round(ecc,2),"    Perihelion=",round(varpi,2)),side=3,font=2)
      abline(h=c(90-obl,obl,0,-1*obl,-90+obl),col="black",lwd=1,lty=4)
      text(50,-90+obl,round(-90+obl,2),font=2)
      text(50,-1*obl,round(-1*obl,2),font=2)
      text(50,0,"0",font=2)
      text(50,obl,round(obl,2),font=2)
      text(50,90-obl,round(90-obl,2),font=2)
     }
     
   if(output & pick==0) return(M) 
   if(output & pick==1) return(data.frame(cbind(long2,MM[,pts])))
   if(output & pick==2) return(data.frame(cbind(lat,MM[pts,])))

### END function insoMap
}

