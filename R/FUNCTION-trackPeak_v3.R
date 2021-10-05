### This function is a component of astrochron: An R Package for Astrochronology
### Copyright (C) 2021 Stephen R. Meyers
###
###########################################################################
### trackPeak: this is a tool to interactively select points to track peak 
###            trajectories on plots, for results from such functions as 
###            eTimeOpt, EHA, eAsm (SRM: December 7, 2017; January 14, 2021;
###                                      August 30, 2021)
###  
###########################################################################

trackPeak <- function (dat,threshold=NULL,pick=T,minVal=NULL,maxVal=NULL,dmin=NULL,dmax=NULL,xmin=NULL,xmax=NULL,ymin=NULL,ymax=NULL,h=6,w=4,ydir=-1,palette=6,ncolors=100,genplot=T,verbose=T)
{
  
  if(verbose) 
    {
      cat("\n---- INTERACTIVELY TRACK POINTS ON PLOT ----\n")
    }
    
# ensure we have a data frame
  dat=data.frame(dat)

# assign sedimentation rates from first column of dat
  sedrates=dat[,1]

  cols=length(dat)
# assign locations for each window (column headers)
  loc=suppressWarnings(as.numeric(substr(colnames(dat[2:cols]),start=2,stop=100)))
# for negative depth/height/time values, "-" has been changed to "."
# this will create NAs. implement modification of fix recommended by Mathieu Martinez
  neg=grepl(".",substr(colnames(dat[2:cols]), start=2,stop=2),fixed=T)
  fixloc=which(neg)
  if(any(neg)) {loc[fixloc]=-1*as.numeric(substr(colnames(dat[(fixloc+1)]),start=3,stop=100))}
# assign r2 values
  sp=as.matrix( dat[2:cols] )

  numrec=length(loc)
  numsed=length(sedrates)
  if(verbose) cat("\n * Number of windows to analyze =",numrec,"\n")
  if(verbose) cat(" * Parameter grid size=",numsed,"\n")

# for analysis
  if(is.null(minVal)) minVal = min(sedrates)
  if(is.null(maxVal)) maxVal = max(sedrates)
  if(is.null(dmin)) dmin = min(loc)
  if(is.null(dmax)) dmax = max(loc)
  
  if(genplot) 
   {

# for plotting
      if(is.null(xmin)) xmin = min(sedrates)
      if(is.null(xmax)) xmax = max(sedrates)
      if(is.null(ymin)) ymin = min(loc)
      if(is.null(ymax)) ymax = max(loc)

# use fields library for access to 'tim.colors', and viridisLite for access to 'viridis'
      if( palette != 1 && palette != 2 && palette != 3 && palette != 4 && palette != 5 && palette != 6) 
        {
          cat("\n**** WARNING: palette option not valid. Will use palette = 6.\n")
          palette = 6
        }

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

# set up device
      dev.new(height=h,width=w)
      par(mfrow=c(1,1))
      xlimset=c(xmin,xmax)

      if (ydir == -1) ylimset=c(ymax,ymin)  
      if (ydir == 1)  ylimset=c(ymin,ymax)
      image(sedrates,loc,sp,xlim=xlimset,ylim=ylimset,col = colPalette,xlab="Parameter",ylab="Depth/Height/Time",main="Click on plot to select peaks")

# end genplot section
   } 

# isolate portion of the results that you want to analyze
   ised= which( (sedrates >= minVal) & (sedrates <= maxVal) )
   sedrates2=sedrates[ised]
   sedrates2=sedrates2[sedrates2 <= maxVal]

   if(verbose) cat("\n * PLEASE WAIT: Sorting windows\n")
 
# perform sorting
  res=rep(NA,numrec)
  for (i in 1:numrec)
    {
# isolate portion of the results that you want to analyze (part 2: 'sp')
      sp2 = sp[ised,i]
      res[i]=peak(cbind(sedrates2,sp2),level=threshold,genplot=F,verbose=F)[2]
    }


# now rearrange results into two column format (this is potentially slow, and should be optimized!)
  outSedrate = 0
  outloc = 0
  for(i in 1:numrec)
   {
# note, is.null sometimes (always?) fails, so use is.na     
     test=(data.frame(res[i]))[1,1]
     if(!is.null(res[i]) || !is.na(test))
      {
        resIn = (data.frame(res[i]))[,1]
        locIn = rep(loc[i],length(resIn))
        outSedrate=append(outSedrate,resIn)
        outloc=append(outloc,locIn)
      } 
   }    

# remove first 'dummy' value
  outloc=outloc[2:length(outloc)]
  outSedrate=outSedrate[2:length(outSedrate)]
  out=data.frame(cbind(outloc,outSedrate))
# Sort to ensure Depth/Height is in increasing order  
#  out <- out[order(out[,1],na.last=NA,decreasing=F),]

## this script modified from '?identify' in R
identifyPch <- function(x, y=NULL, n=length(x), pch=19, ...)
{
    xy <- xy.coords(x, y); x <- xy$x; y <- xy$y
    sel <- rep(FALSE, length(x)); res <- integer(0)
    while(sum(sel) < n) {
        ans <- identify(x[!sel], y[!sel], n=1, plot=F, ...)
        if(!length(ans)) break
        ans <- which(!sel)[ans]
        points(x[ans], y[ans], pch = pch, col="white")
        sel[ans] <- TRUE
        res <- c(res, ans)
    }
    res
}

  if(is.null(threshold) && genplot && verbose) 
    {
      cat("\n**** WARNING: By default all peak maxima are reported and plotted.\n") 
      cat("              If there are many peaks, this will cause your plot to be mostly white.\n")
      cat("              Set threshold to cull peaks.\n") 
    }  


# Now plot, identifying location of peaks
  if(genplot)
    {
       par(new=T)
       plot(out[,2],out[,1],xlim=xlimset,ylim=ylimset,xaxs="i",yaxs="i",yaxt='n',bty='n',ylab="",xlab="",pch=1,col="white")
# if interactive point identification selected
       if(pick)
        {
          cat("\n---- INTERACTIVELY SELECT POINTS ----\n")
          cat("\ *****  Select path by clicking  *****\n")
          cat(" Stop by pressing ESC-key (Mac) or STOP button (Windows)\n")
          pts <- identifyPch(out[,2], out[,1])
          out=data.frame(cbind(out[pts,1],out[pts,2]))
        }  
# end genplot section
    }

  cat("\n * Peaks identified.\n")

# sort out to ensure increasing depth/height
  out <- out[order(out[,1],na.last=NA,decreasing=F),]
  colnames(out)[1] = 'Depth/Height/Time'
  colnames(out)[2] = 'Parameter'

  return(out)

### END function trackPeak
}
