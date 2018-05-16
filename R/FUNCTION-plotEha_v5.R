### This function is a component of astrochron: An R Package for Astrochronology
### Copyright (C) 2018 Stephen R. Meyers
###
###########################################################################
### plotEha: create time-frequency plots from eha results using image 
###          or image.plot (SRM: Aug. 8, 2013; April 20, 2014; January 16, 2015
###                              February 10, 2017; January 30, 2018)
###                                                                 
###########################################################################

plotEha <- function (spec,xmin=NULL,xmax=NULL,ymin=NULL,ymax=NULL,h=6,w=4,ydir=1,pl=0,norm=NULL,palette=1,centerZero=T,ncolors=100,colorscale=F,xlab=NULL,ylab=NULL,filetype=0,output=T,verbose=T)
{  

# uses fields library for tim.colors

  if(verbose) cat("\n---- PLOTTING EHA OUTPUT ----\n")

# ensure we have a data frame
  spec=data.frame(spec)

# error checking
  if( palette != 1 && palette != 2 && palette != 3 && palette != 4 && palette != 5) 
     {
       cat("\n**** WARNING: palette option not valid. Will use palette = 1.\n")
       palette = 1
     }
     
  if (is.null(xlab)) xlab = c("Frequency")
  if (is.null(ylab)) ylab = c("Location")

# assign frequencies from first column of spec
  freq=spec[,1]

  cols=length(spec)
# assign locations for each spectrum (column headers)
  loc=suppressWarnings(as.numeric(substr(colnames(spec[2:cols]),start=2,stop=100)))
# for negative depth/height/time values, "-" has been changed to "."
# this will create NAs. implement modification of fix recommended by Mathieu Martinez
  neg=grepl(".",substr(colnames(spec[2:cols]), start=2,stop=2),fixed=T)
  fixloc=which(neg)
  if(any(neg)) {loc[fixloc]=-1*as.numeric(substr(colnames(spec[(fixloc+1)]),start=3,stop=100))}
# assign specta (amplitude, power, or probability)
  sp=as.matrix( spec[2:cols] )

  numrec=length(loc)
  numfreq=length(freq)
  if(verbose) cat("\n * Number of spectra =",numrec,"\n")
  if(verbose) cat(" * Number of frequencies per spectrum =",numfreq,"\n")

  if(is.null(norm) && pl == 3) stop("**** ERROR: variable norm is not defined! Terminating now!")
  if(!is.null(norm)) norm=data.frame(norm)

# use fields library for access to 'tim.colors' only

# for plotting
  if(is.null(xmin)) xmin = min(freq)
  if(is.null(xmax)) xmax = max(freq)
  if(is.null(ymin)) ymin = min(loc)
  if(is.null(ymax)) ymax = max(loc)
  if(pl==0) sp_plot=sp
  if(pl==1) sp_plot=log(sp)
  if(pl==2) 
       {
           normAmp = apply(sp,2,max)
           sp_plot=t(sp)/normAmp
           sp_plot=t(sp_plot)
       }
  if(pl==3)
       {
           sp_plot=t(sp)/norm[,1]
           sp_plot=t(sp_plot)
       }           

# remove -Inf values, if present
  if(min(sp_plot== -Inf)) 
    {
      if(verbose) cat("\n**** WARNING: -Inf values will be replaced with",-1*.Machine$double.xmax,", for plotting.\n")
      sp_plot[sp_plot==-Inf] <- -1*.Machine$double.xmax
    }
    
# set color palette
#  rainbow colors
 if(palette == 1) colPalette = tim.colors(ncolors)
#  grayscale
 if(palette == 2) colPalette = gray.colors(n=ncolors,start=1,end=0,gamma=1.75)
#  dark blue scale (from larry.colors)
# if(palette == 3) colPalette = colorRampPalette(c("white","#023858"))(ncolors)
 if(palette == 3) colPalette = colorRampPalette(c("white","royalblue"))(ncolors)
#  red scale
 if(palette == 4) colPalette = colorRampPalette(c("white","red2"))(ncolors)
#  blue to red plot
# if(palette == 5) colPalette = append(colorRampPalette(c("#023858","white"))(ncolors/2),colorRampPalette(c("white","red2"))(ncolors/2))
 if(palette == 5) colPalette = append(colorRampPalette(c("royalblue","white"))(ncolors/2),colorRampPalette(c("white","red2"))(ncolors/2))

# this centers the color scale on zero (an equal number of color divisions above and below zero), as recommended by Nature (Feb. 10, 2017)
 if(max(sp_plot) > 0 && min(sp_plot) < 0 && centerZero)
  {
    break1=seq(min(sp_plot),0,length.out=(ncolors/2)+1)
    break2=seq(min(sp_plot[sp_plot>0]),max(sp_plot),length.out=(ncolors/2))
    breaksPwr=append(break1,break2)
  }           
if(!centerZero || max(sp_plot) <= 0 || min(sp_plot) >= 0) breaksPwr=seq(min(sp_plot),max(sp_plot),length.out=(ncolors)+1)  

# set up device
  if(filetype==0) dev.new(height=h,width=w)
  if(filetype==1) pdf("plotEha.pdf",height=h,width=w)
  if(filetype==2) jpeg("plotEha.jpeg",height=h,width=w,units="in",res=250)
  if(filetype==3) png(filename = "plotEha.png",width=w,height=h,units="in",res=250)
  if(filetype==4) tiff(filename = "plotEha.tiff",width=w,height=h,units="in",res=250)


  par(mfrow=c(1,1))
  xlimset=c(xmin,xmax)

  if (ydir == -1) 
    {
# in this case, reset ylim range.
# note that useRaster=T is not a viable option, as it will plot the results backwards, even though the
#  y-axis scale has been reversed!  This option will result in a slower plotting time.
        ylimset=c(ymax,ymin)
        if(colorscale) image.plot(freq,loc,sp_plot,xlim=xlimset,ylim=ylimset,col=colPalette,breaks=breaksPwr,xlab=xlab,ylab=ylab)
        if(!colorscale) image(freq,loc,sp_plot,xlim=xlimset,ylim=ylimset,col=colPalette,breaks=breaksPwr,xlab=xlab,ylab=ylab)
     }

   if (ydir == 1) 
     {
# useRaster=T results in a faster plotting time.
        ylimset=c(ymin,ymax)
        if(colorscale) image.plot(freq,loc,sp_plot,xlim=xlimset,ylim=ylimset,col=colPalette,breaks=breaksPwr,useRaster=T,xlab=xlab,ylab=ylab)       
        if(!colorscale) image(freq,loc,sp_plot,xlim=xlimset,ylim=ylimset,col=colPalette,breaks=breaksPwr,useRaster=T,xlab=xlab,ylab=ylab)
      }

     if(filetype != 0) dev.off()
     if(output && pl == 2) return(data.frame(normAmp))

### END function plotEha
}