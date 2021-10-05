### This function is a component of astrochron: An R Package for Astrochronology
### Copyright (C) 2021 Stephen R. Meyers
###
###########################################################################
### function mwinGrid - (SRM: August 29, 2020; January 14, 2021)
###
### determine start and end points for a moving window of fixed
###  duration (e.g., 500 kiloyears), using a fixed sampling grid (space or time)
###  that is evenly spaced.
##  for use with evenly or unevenly sampled data.
### no duplcate values or empty entries (NA) permitted.
###########################################################################

mwinGrid <- function (dat,win,step,start=NULL,end=NULL,verbose=T)
{

   if(verbose) cat("\n----- DETERMINING PARAMETERS OF MOVING WINDOW -----\n")

# make sure the input series is a data frame
   dat=data.frame(dat)

   ipts <- length(dat[,1]) 
   if(verbose) cat(" * Number of data points=", ipts,"\n")

# sort to ensure increasing depth/height/time
   if(verbose) cat(" * Sorting\n")
   dat <- dat[order(dat[,1],na.last=NA,decreasing=F),]

# error checking
   if(length(dat[,1]) != ipts) 
     {
       cat("\n**** ERROR: empty entires found\n")
       stop("**** TERMINATING NOW!")
     }

   dtest <- dat[2:ipts,1]-dat[1:(ipts-1),1]
   if(any(dtest == 0)) 
     {
       cat("\n**** ERROR: duplicate stratigraphic levels found.\n")
       stop("**** TERMINATING NOW!")
     }

   x <- dat[,1]

# middle of first window
   if(is.null(start)) start = x[1]
   begin = start + win/2
# middle of the last window
   if(is.null(end)) end = x[ipts]
   finish = end - win/2
   outGrid=seq(from=begin,to=finish,by=step)
   iGridPts=length(outGrid)


# from FORTRAN: subroutine mwinGrid_R(ipts,x,start,iGridPts,step,winsize,n1,n2)
getwin <- function (ipts,x,start,iGridPts,step,win)
 {
    F_dat = .Fortran( 'mwingrid_r',
    
    ipts=as.integer(ipts),x=as.double(x),start=as.double(start),
    iGridPts=as.integer(iGridPts),step=as.double(step),winsize=as.double(win),
    
    n1=integer(iGridPts),n2=integer(iGridPts)
    )

# return the results
    return(F_dat)
 }

   out <- getwin(ipts,x,start,iGridPts,step,win)

   out2=data.frame(cbind(out$n1,out$n2,outGrid))
   colnames(out2) = c("n1","n2","center")
   
   
   if(verbose)
    {
      if(any(out2[,1]==out2[,2])) cat("\n**** WARNING: some windows contain only one data point!\n")
      if(any(out2[,1]>out2[,2])) cat("\n**** WARNING: some windows contain no data points!\n\n")
    }
       
   if(any(out2[,1]>out2[,2])) 
    {
      makeNA=which(out2[,1]>out2[,2])
      out2[makeNA,1] <- NA
      out2[makeNA,2] <- NA
    }   
       
       
   return(out2)

### END function mwinGrid
}
