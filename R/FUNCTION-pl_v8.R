### This function is a component of astrochron: An R Package for Astrochronology
### Copyright (C) 2020 Stephen R. Meyers
###
###########################################################################
### pl: set up plotting parameters (SRM: March 29, 2012; June 27, 2013; 
###                                      July 31, 2013; September 16, 2013
###                                      January 13, 2014; January 14, 2014;
###                                      February 13, 2014; June 30, 2014;
###                                      February 3, 2015; December 5, 2020;
###                                      December 14, 2022)
###
###########################################################################

pl <- function (n=NULL,r=NULL,c=NULL,h=NULL,w=NULL,mar=NULL,file=NULL,title=NULL)
{

if(is.null(h) || is.null(w)) 
 {
   cat("\n**** WARNING: using default height and width of 7 inches.\n")
   h=7
   w=7
 }

if(is.null(mar)) mar=c(5.1,4.1,4.1,2.1)
if(!is.null(mar) && length(mar)!=4) 
 {
   cat("\n**** WARNING: mar not specified correctly.  Will use default values.\n")
   mar=c(5.1,4.1,4.1,2.1)
 }  

if(is.null(file)) 
 {
# use default window
   dev.new(height=h,width=w,title=title)
  } 
if(!is.null(file)) 
 {
   if(grepl(".pdf",file,fixed=TRUE)) pdf(file,height=h,width=w)
   if(grepl(".jpg",file)) jpeg(file,height=h,width=w,units="in",res=300)
   if(grepl(".jpeg",file)) jpeg(file,height=h,width=w,units="in",res=300)
   if(grepl(".png",file)) png(file,height=h,width=w,units="in",res=300)
   if(grepl(".tiff",file)) tiff(file,height=h,width=w,units="in",res=300)
   if(grepl(".bmp",file)) bmp(file,height=h,width=w,units="in",res=300)
 }
 
if(is.null(n))  
 { 
  if(is.null(r) || is.null(c))
   {
# by default, set up for one plot
     if(is.null(r) && is.null(c))
      {
        cat(" * No parameters specified.  Will set up device for one plot.\n")
        r = 1
        c = 1
       }
     else
       {
         r=max(r,c)
         c=max(r,c)
         cat(" * Missing input parameter. Will set up device for",r,"rows,",c,"columns.\n")
       }
   }
  par(mfrow=c(r,c))
 }

# if n is specified, it takes precedence  
if(!is.null(n))
 {
   if(!is.null(r)) cat("**** WARNING: n takes precedence, r will be ignored.\n")
   if(!is.null(c)) cat("**** WARNING: n takes precedence, c will be ignored.\n")
   if(n >25) 
    {
      cat("**** WARNING: Maximum n=25")
      n = 25
    }
   if(n == 1) par(mfrow=c(1,1))
   if(n == 2) par(mfrow=c(2,1))
   if(n == 3) par(mfrow=c(3,1))
   if(n == 4) par(mfrow=c(2,2))
   if(n == 5) par(mfrow=c(3,2))
   if(n == 6) par(mfrow=c(3,2))
   if(n == 7 || n == 8) par(mfrow=c(4,2))
   if(n == 9) par(mfrow=c(3,3))
   if(n == 10 || n==12) par(mfrow=c(4,3))
   if(n > 12 && n<17) par(mfrow=c(4,4))
   if(n > 16 && n<21) par(mfrow=c(5,4))
   if(n > 20 && n<26) par(mfrow=c(5,5))
 }

if(!is.null(file)) cat(" * Be sure to close device after plotting is complete: dev.off()")

    par(mar=mar)

### END function pl
}
