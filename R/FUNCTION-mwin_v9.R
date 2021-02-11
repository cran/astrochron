### This function is a component of astrochron: An R Package for Astrochronology
### Copyright (C) 2021 Stephen R. Meyers
###
###########################################################################
### function mwin - (SRM: March 15, 2013; May 20, 2013; June 5, 2013; 
###                       May 24, 2017; May 31, 2017; June 12-14, 2017;
###                       January 14, 2021)
###
### determine start and end points for a moving window of fixed
### duration (e.g., 500 kiloyears). for use with unevenly sampled data.
### no duplcate values or empty entries (NA) permitted.
###########################################################################

mwin <- function (dat,win,conv=1,verbose=T)
{

   if(verbose) cat("\n----- DETERMINING PARAMETERS OF DYNAMIC MOVING WINDOW -----\n")

# make sure the input series is a data frame
   dat=data.frame(dat)

   ipts <- length(dat[,1]) 
   if(verbose) cat(" * Number of data points=", ipts,"\n")

# error checking
   if(conv != 1 && conv != 2 && conv != 3) 
     {
       cat("\n**** ERROR: conv must be set to 1, 2, or 3\n")
       stop("**** TERMINATING NOW!")
     }

# sort to ensure increasing depth/height/time
   if(conv == 1 || conv == 2)
    {
     if(verbose) cat(" * Sorting\n")
     dat <- dat[order(dat[,1],na.last=NA,decreasing=F),]
    }

# flip stratigraphic series around for getwin function
   if(conv == 3)
    {
     if(verbose) cat(" * Sorting\n")
     dat <- dat[order(dat[,1],na.last=NA,decreasing=T),]
     dat[1]=-1*dat[1]
    }

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

getwin1 <- function (ipts,win,x)
 {
    F_dat = .Fortran( 'mwincenter_r',
    
    ipts=as.integer(ipts),winsize=as.double(win),x=as.double(x),
    
    npts=integer(1),n1=integer(ipts),n2=integer(ipts),avex=double(ipts),
    midx1=double(ipts),midx2=double(ipts) 
    )

# return the results
    return(F_dat)
 }

getwin2 <- function (ipts,win,x)
 {
    F_dat = .Fortran( 'mwin_r',
    
    ipts=as.integer(ipts),winsize=as.double(win),x=as.double(x),
    
    npts=integer(1),n1=integer(ipts),n2=integer(ipts),avex=double(ipts),
    midx1=double(ipts),midx2=double(ipts) 
    )

# return the results
    return(F_dat)
 }

# call the function
   if(conv == 1) out <- getwin1(ipts,win,x)
   if(conv == 2 || conv == 3) out <- getwin2(ipts,win,x)

   npts=out$npts
   out2=data.frame(cbind(out$n1[1:npts],out$n2[1:npts],out$avex[1:npts],out$midx1[1:npts],out$midx2[1:npts]))
   colnames(out2) = c("n1","n2","avex","center","midpoint")
   
   if(conv == 3)
    {
     out2[3] = -1*out2[3]
     out2[4] = -1*out2[4]
     out2[5] = -1*out2[5]
# correct the indices so they correspond to the original input data
# note that n1 and n2 must be swapped
     out2[,1] = ipts - out$n2[1:npts] + 1
     out2[,2] = ipts - out$n1[1:npts] + 1      
# resort
     out2 <- out2[order(out2[,4],na.last=NA,decreasing=F),]    
    } 

   if(verbose)
    {
      if(any(out2[,1]==out2[,2])) cat("\n**** WARNING: some windows contain only one data point!\n")
    }
       
return(out2)

### END function mwin
}
