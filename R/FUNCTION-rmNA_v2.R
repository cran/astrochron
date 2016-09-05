### This function is a component of astrochron: An R Package for Astrochronology
### Copyright (C) 2016 Stephen R. Meyers
###
###########################################################################
### rmNA: remove NA's- (SRM: May 11, 2016; July 22, 2016)
###
###########################################################################

rmNA <- function (dat, genplot=T, verbose=T)
{

  if(verbose) cat("\n----- REMOVING STRATIGRAPHIC LEVELS THAT CONTAIN 'NA' ENTRIES -----\n")

  dat <- data.frame(dat)
  ipts <- length(dat[,1]) 
  if(verbose) cat(" * Number of rows (stratigraphic levels)=", ipts,"\n")
  cols=length(dat)
  if(verbose) cat(" * Number of columns=", cols,"\n")
 
  numNA=sum(is.na(dat))
  if(verbose) cat(" * Number of NA values identified=", numNA,"\n")
  
# processes each stratigraphic level (row), delete if it has an NA in any column
  dat=na.omit(dat)
  
  newpts=length(dat[,1])
  if(verbose) cat(" * Number of rows (stratigraphic levels) following culling=",newpts,"\n")
   
 if(cols == 2 && genplot)
  {
### plots
   par(mfrow=c(2,2))
   plot(dat,cex=0.5,xlab="Location",ylab="Value",main="Data Series"); lines(dat)
### plot the denisty and the histogram together
   hist(dat[,2],freq=F,xlab="Value",main="Distribution of isoloated values"); lines(density(dat[,2], bw="nrd0"),col="red")
### boxplot
   boxplot(dat[,2],ylab="Value",main="Boxplot for isolated values")
### Normal probabilty plot (Normal Q-Q Plot)
   qqnorm(dat[,2]); qqline(dat[,2], col="red")
  }
  
  return(dat)

### END function rmNA
}
