### This function is a component of astrochron: An R Package for Astrochronology
### Copyright (C) 2021 Stephen R. Meyers
###
###########################################################################
### tracePlot: this is a tool to interactively trace peak trajectories on 
###            plots, for results from such functions as eTimeOpt, eha, eAsm
###            (SRM: December 7, 2017; January 14, 2021; August 30, 2021)
###  
###########################################################################

tracePeak <- function (dat,color=2,h=6,w=4,xmin=NULL,xmax=NULL,ymin=NULL,ymax=NULL,ydir=-1,palette=6,ncolors=100,path=1)
{

  cat("---- INTERACTIVELY TRACE TRAJECTORY ON PLOT ----\n")
  cat("\n  *****  Select path by clicking     *****\n")
  cat("  Stop by pressing ESC-key (Mac) or STOP button (Windows)\n")

# use fields library for access to 'tim.colors', and viridisLite for access to 'viridis'
  if( palette != 1 && palette != 2 && palette != 3 && palette != 4 && palette != 5 && palette != 6) 
    {
       cat("\n**** WARNING: palette option not valid. Will use palette = 6.\n")
       palette = 6
    }

  dat=data.frame(dat)

# assign sedimentation rates from first column of dat
  sedrates=dat[,1]
  rows=length(sedrates)

  cols=length(dat)
# assign locations for each window (column headers)
  loc=suppressWarnings(as.numeric(substr(colnames(dat[2:cols]),start=2,stop=100)))
# for negative depth/height/time values, "-" has been changed to "."
# this will create NAs. implement modification of fix recommended by Mathieu Martinez
  neg=grepl(".",substr(colnames(dat[2:cols]), start=2,stop=2),fixed=T)
  fixloc=which(neg)
  if(any(neg)) {loc[fixloc]=-1*as.numeric(substr(colnames(dat[(fixloc+1)]),start=3,stop=100))}
# assign specta
  sp=as.matrix( dat[2:cols] )

  if(is.null(xmin)) xmin = min(sedrates)
  if(is.null(xmax)) xmax = max(sedrates)
  if(is.null(ymin)) ymin = min(loc)
  if(is.null(ymax)) ymax = max(loc)
  if(path==1) pltype="o"
  if(path==2) pltype="l"
  if(path==3) pltype="p"

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

# set up plot
  dev.new(height=h,width=w)
  par(mfrow=c(1,1))
  xlimset=c(xmin,xmax)

  if (ydir == -1) ylimset=c(ymax,ymin)
  if (ydir == 1) ylimset=c(ymin,ymax)
  image(sedrates,loc,sp,xlim=xlimset,ylim=ylimset,col = colPalette,xlab="Parameter",ylab="Depth/Height/Time",main="Click on plot to define trajectory")       
     
# Now overlay x-y plot for graphical interface
  par(new=T)
  plot(-1,-1,xlim=xlimset,ylim=ylimset,xaxs="i",yaxs="i",yaxt='n',bty='n',ylab="",xlab="")
# transparent black
  if (color == 1) setcolor="#00000070"
# transparent white: alpha (transparency) of 210 in hexadecimal = D2       
  if (color == 2) setcolor="#FFFFFFD2"
# transparent yellow
  if (color == 3) setcolor="#FFFF00D2"
  ff = locator(n = length(loc), type = pltype ,col=setcolor,lwd=3)

  out = data.frame(cbind(ff$y,ff$x))
  colnames(out)[1] = 'Depth/Height/Time'
  colnames(out)[2] = 'Parameter'
       
# Sort to ensure Depth/Height is in increasing order
  out <- out[order(out[,1],na.last=NA,decreasing=F),]

  cat("\n * Tracing complete.\n")

  return(out)

### END function tracePeak
}
