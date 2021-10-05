### This function is a component of astrochron: An R Package for Astrochronology
### Copyright (C) 2021 Stephen R. Meyers
###
###########################################################################
### mwStatsGrid : moving window average, median, and standard deviation, 
###               allowing for dynamic adustment of window so that it has a 
###               constant duration in time or space. This version will conduct
###               the analysis on an evenly spaced spatial/temporal grid. 
###                 (SRM: August 31, 2020; September 20, 2020; January 14, 2021;
###                       August 30, 2021) 
###
###########################################################################
# modified from FUNCTION-mwStats_v6.R
# may want to add check for the case when no points are in the window, and option to
#  change direction ('ydir' option) when plotting kernel density estimates.

mwStatsGrid <- function (dat,cols=NULL,win=NULL,step=NULL,start=NULL,end=NULL,output=T,norm=F,palette=6,ncolors=100,genplot=1,verbose=T)

{

  if(verbose) cat("\n----- CALCULATING MOVING AVERAGE, MEDIAN AND VARIANCE USING DYNAMIC WINDOW-----\n")

# ensure we have a data frame
  dat=data.frame(dat)
  
  if(is.null(cols) && length(dat) != 2)
    {
      if(verbose) cat("\n**** ERROR: number of variables > 1.  Please specify which variable to analyze.\n")  
      stop("**** TERMINATING NOW")
    }
# no need to do anything if is.null(cols) && length(dat) == 2
  if(!is.null(cols) && length(dat) == 2) 
    {
     if(verbose) cat("\n**** WARNING: You only have two columns, will ignore cols.\n") 
    } 
  if(!is.null(cols) && length(cols) != 1 && length(dat) > 2) 
    {
      if(verbose) cat("\n**** WARNING: cols not specified correctly, will use data from column 2.\n")  
# reassign dat
      dat <- data.frame( cbind(dat[,1],dat[,2]) )
    }  
  if(!is.null(cols) && length(cols) == 1 && length(dat) > 2) 
    {
# reassign dat
      dat <- data.frame( cbind(dat[,1],dat[,cols[1]]) )
    }  

### sort to ensure increasing depth/height/time
   if(verbose) cat(" * Sorting into increasing order, removing empty entries\n")
   dat <- dat[order(dat[,1],na.last=NA,decreasing=F),]

   ipts <- length(dat[,1]) 
   if(verbose) cat(" * Number of data points=", ipts,"\n")  

   if(is.null(win)) 
    {
      win=(dat[ipts,1]-dat[1,1])/10
      if(verbose) cat(" * Setting window size to default value of", win,"\n")  
    }  
  
   if(win >= (dat[ipts,1]-dat[1,1])) 
    {
      if(verbose) cat("\n**** ERROR: window length is >= data series duration.\n")  
      stop("**** TERMINATING NOW")    
    } 

   if(is.null(step)) 
    {
      step=(dat[ipts,1]-dat[1,1])/100
      if(verbose) cat(" * Setting window step size to default value of", step,"\n")  
    }  

   
# determine location of each window
  loc <- mwinGrid(dat,win,step,start,end,verbose=F)
  
  npts=length(loc$n1)
  res=double(npts)
  res2=double(npts)
  res3=double(npts)
  n=double(npts)
  
# now loop over all windows
  for (i in 1:npts)
   {
# shuffle in data
     if(loc$n1[i]==loc$n2[i]) 
      {
       cat("**** WARNING: Only one data point in window",i,":",loc$center[i],"\n")
        y <- dat[loc$n1[i],2]
        n[i]=1
        res[i] = y
        res2[i] = y
        res3[i] = var(y)
      }
     
     if(loc$n1[i]!=loc$n2[i])
      {   
        y <- dat[loc$n1[i]:loc$n2[i],2]
# number of data points in window
        n[i]=loc$n2[i]-loc$n1[i]+1
        res[i] = mean(y)
        res2[i] = median(y)
        res3[i] = var(y)
      }  
    }

 if(genplot==1)
   {
     dev.new(height = 7.7, width = 5.7, title = "mwStatsGrid Results")
     par(mfrow=c(4,1),mar=c(3.1, 4.1, 4.1, 2.1))
     
     ymin=min(res,res2)
     ymax=max(res,res2)
     xmin=min(dat[,1])
     xmax=max(dat[,1])
     plot(dat[,1],dat[,2],type="l",xlab="",ylab="",main="Stratigraphic Series",bty="n",xlim=c(xmin,xmax))
     mtext("Location",side=1,line=2,cex=0.7)
     mtext("Value",side=2,line=2,cex=0.7)

     par(mar=c(4.1, 4.1, 3.1, 2.1))
     plot(loc$center,res,ylim=c(ymin,ymax),xlim=c(xmin,xmax),type="l",lwd=2,xlab="",ylab="",main="Moving Window Average (black) and Median (red)",bty="n")
     lines(loc$center,res2,col="red",lwd=2)

     mtext("Center of window",side=1,line=2,cex=0.7)
     mtext("Value",side=2,line=2,cex=0.7)

     plot(loc$center,res3,xlim=c(xmin,xmax),type="l",xlab="",ylab="",main="Variance",lwd=2,bty="n")
     mtext("Center of window",side=1,line=2,cex=0.7)
     mtext("Variance",side=2,line=2,cex=0.7)
     
     plot(loc$center,n,xlim=c(xmin,xmax),type="l",xlab="",ylab="",main="Number of Data Points in Window",bty="n")
     mtext("Center of window",side=1,line=2,cex=0.7)
     mtext("# Points",side=2,line=2,cex=0.7)
   }
   
  if(genplot >= 2)
     {

       if( palette != 1 && palette != 2 && palette != 3 && palette != 4 && palette != 5 && palette != 6) 
         {
           cat("\n**** WARNING: palette option not valid. Will use palette = 6.\n")
            palette = 6
         }

# uses fields library for plotting
# get common grid for all density results
       gridDen=density(dat[,2])$x
       nGridDen=length(gridDen)       
# set up matrices for results (gridDen already defined)
       height<-loc$center
       densityOut<-double(nGridDen*npts)
       dim(densityOut)<-c(nGridDen,npts)

# now loop over all windows
       for (i in 1:npts)
          {
            if(loc$n1[i]!=loc$n2[i])
              {   
                y <- dat[loc$n1[i]:loc$n2[i],2]
# number of data points in window
                densityOut[,i] = density(y,bw="nrd0",na.rm=T,from=gridDen[1],to=gridDen[nGridDen],n=nGridDen)$y
                if(norm) densityOut[,i] = densityOut[,i]/max(densityOut[,i])
              }  
          }  
# plot
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

       dev.new(title=paste("mwStatsGrid results"),height=6,width=4.5)
       xlimset=c(min(gridDen),max(gridDen))
       ylimset=c(min(loc$center),max(loc$center))
       image.plot(gridDen,height,densityOut,col=colPalette,xlim=xlimset,ylim=ylimset,xlab="Value",ylab="Position",main="Moving Window Kernel Density Plots")

       if(genplot >=3)
         {
           par(new=T)
           if(genplot == 3) 
             {
               plot(res2,loc$center,type="l",xlim=xlimset,ylim=ylimset,xaxs="i",xaxt="n",yaxs="i",yaxt="n",bty="n",ylab="",xlab="",lwd=3,col="white")
               lines(res2,loc$center,lwd=2,col="black")
             }  
           if(genplot == 4) 
             {
               plot(res,loc$center,type="l",xlim=xlimset,ylim=ylimset,xaxs="i",xaxt="n",yaxs="i",yaxt="n",bty="n",ylab="",xlab="",lwd=3,col="white")
               lines(res,loc$center,lwd=2,col="black")
             }  
         }
     }

    if(output)
     {
       out = data.frame (cbind (loc$center,res,res2,res3,n) ) 
       colnames(out)[1] <- 'Center_win'
       colnames(out)[2] <- 'Average'
       colnames(out)[3] <- 'Median'
       colnames(out)[4] <- 'Variance'
       colnames(out)[5] <- 'Points'
       return(out)
     }  

### END function mwStatsGrid
}