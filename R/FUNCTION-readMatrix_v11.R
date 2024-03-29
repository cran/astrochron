### This function is a component of astrochron: An R Package for Astrochronology
### Copyright (C) 2023 Stephen R. Meyers
###
###########################################################################
### readMatrix - (SRM: May 2, 2014; September 28, 2014; January 22, 2015; 
###               March 1, 2015; August 6, 2015; September 10, 2015; 
###               September 16, 2015; November 18, 2016; July 26, 2017;
###               December 15, 2017; December 5, 2020; July 27, 2023)
###
### Read a matrix from file 
### The file can include any number of columns, each representing a different variable.
###########################################################################

readMatrix <- function (file=NULL,d=1,h="auto",skip=0,output=1,check=T,genplot=F,verbose=T)
{
 if(verbose)
  {
   cat("\n----- READ MATRIX FROM DATA FILE -----\n")
   cat("\nThe following options are selected:\n")
   if(d==0) cat(" * What type of column delimiter are you using?: Tab\n")
   if(d==1) cat(" * What type of column delimiter are you using?: Comma\n")
   if(d==2) cat(" * What type of column delimiter are you using?: Semicolon\n")
   if(h=="yes") cat(" * Does your file have column titles/headers?: yes\n")
   if(h=="no") cat(" * Does your file have column titles/headers?: no\n")
   if(h=="auto") cat(" * Does your file have column titles/headers?: auto detect\n")
   if(skip>0) cat(" * Skipping", skip, "line(s) of file\n")
      
# pause so there is time to ouput text to screen.
   Sys.sleep(0.5)
  }

### if file name and path set
   if(!is.null(file)) filen <- file
### if not, get file interactively
   if(is.null(file))
    {
      cat("\n  <PLEASE CHOOSE YOUR FILE>\n")
      filen <- file.choose()
    }  

  if(h=="no")
     {
       if(d == 0) {dat <- read.table (filen,header=F,skip=skip)}
       if(d == 1) {dat <- read.table (filen,header=F,sep=",",skip=skip)}
       if(d == 2) {dat <- read.table (filen,header=F,sep=";",skip=skip)}
       xlab="Value"
       ylab="Value"
     }

  if(h=="yes")
     {
       if(d == 0) {dat <- read.table (filen,header=T,skip=skip)}
       if(d == 1) {dat <- read.table (filen,header=T,sep=",",skip=skip)}
       if(d == 2) {dat <- read.table (filen,header=T,sep=";",skip=skip)}
       xlab=names(dat[1])
       if (length(dat) > 1) ylab=names(dat[2])
     }

# auto detect column titles, for tab-delimited files
   if(h=="auto" && d == 0) 
      {
        dat <- read.table (filen,nrows=1,colClasses="character",skip=skip)
        titles=is.na(suppressWarnings(as.numeric(dat[1,1])))
# no column titles
        if (!titles )  
         {
           if(verbose) cat("\n * No column titles/headers detected\n")
           dat <- read.table (filen,header=F,skip=skip)
           xlab="Value"
           if (length(dat) > 1) ylab="Value"
         }
# column titles        
        if (titles )  
         {
           if(verbose) cat("\n * Column titles/headers detected\n")
           dat <- read.table (filen,header=T,skip=skip)
           xlab=names(dat[1])
           if (length(dat) > 1) ylab=names(dat[2]) 
         }
      }

# auto detect column titles, for csv files
   if(h=="auto" && d == 1) 
      {
        dat <- read.table (filen,nrows=1,colClasses="character",sep=",",skip=skip)
        titles=is.na(suppressWarnings(as.numeric(dat[1,1])))
# no column titles
        if (!titles )  
         {
           if(verbose) cat("\n * No column titles/headers detected\n")
           dat <- read.table (filen,header=F,sep=",",skip=skip)
           xlab="Location"
           if (length(dat) > 1) ylab="Value"         
         }
# column titles        
        if (titles)  
         {
          if(verbose) cat("\n * Column titles/headers detected\n")
          dat <- read.table (filen,header=T,sep=",",skip=skip)
          xlab=names(dat[1])
          if (length(dat) > 1) ylab=names(dat[2])    
         }
      }

# auto detect column titles, for semicolon-delimited files
   if(h=="auto" && d == 2) 
      {
        dat <- read.table (filen,nrows=1,colClasses="character",sep=";",skip=skip)
        titles=is.na(suppressWarnings(as.numeric(dat[1,1])))
# no column titles
        if (!titles )  
         {
           if(verbose) cat("\n * No column titles/headers detected\n")
           dat <- read.table (filen,header=F,sep=";",skip=skip)
           xlab="Location"
           if (length(dat) > 1) ylab="Value"         
         }
# column titles        
        if (titles)  
         {
          if(verbose) cat("\n * Column titles/headers detected\n")
          dat <- read.table (filen,header=T,sep=";",skip=skip)
          xlab=names(dat[1])
          if (length(dat) > 1) ylab=names(dat[2])    
         }
      }
     
   npts <- length(dat[,1]) 
   if(verbose) cat("\n * Number of rows=", npts,"\n")

   cols <- length(dat[1,])
   if(verbose) cat(" * Number of columns=",cols,"\n")

  if(check)
   {   
# check to see if any of the rows are all NA entries
#  logical will default to FALSE
   delRow <- logical(npts)
   for (i in 1:npts) delRow[i]=all(is.na(dat[i,]))
# delete rows that have all NA entries  
   if(any(delRow))
    { 
      dat <- dat[!delRow,]
      npts <- length(dat[,1]) 
      if(verbose) 
       {
         cat("\n * Some rows contain all NA entries, and will be removed\n")
         cat(" * New number of rows=", npts,"\n")
       }
   }

# check to see if any of the columns are all NA entries
   delCol<-logical(cols)
   for (i in 1:cols) delCol[i]=all(is.na(dat[,i]))
# delete rows that have all NA entries  
   if(any(delCol))
    { 
      dat<-dat[,!delCol]
      cols <- length(dat[1,]) 
      if(verbose) 
       {
         cat("\n * Some columns contain all NA entries, and will be removed\n")
         cat(" * New number of columns=", cols,"\n")
       }
   }

# determine if any missing entries (NA) still exist (this may be the case when cols>2)
     numNA=sum(is.na(dat))
     if(numNA > 0 && verbose) 
       {
         cat("\n  WARNING:", numNA,"empty (NA) entries are still present in your data frame.\n")
       } 
     }

if(genplot) 
  {
    if(cols >1) nrow = ceiling(sqrt(cols-1))
    if(cols == 1) nrow = 1
    ncol = nrow
    par(mfrow = c(nrow, ncol))
    for (i in 1:cols) 
      {
        xlab = names(dat)[i]
        plotdat = subset(dat[,i], !(dat[,i] == "NA"))
        if( inherits(dat[,i],"numeric") )
          {
            hist(plotdat,freq=F,xlab=xlab,main=""); lines(density(plotdat, bw="nrd0"),col="red")
          } 
        if( !inherits(dat[,i],"numeric") )
          {
        if(verbose) cat("\n  WARNING: column",i,"contains non-numeric values.  Will NOT generate a plot.\n")
          }
      }
  }
    
# output headers
   if(verbose)
    {
      cat("\n * First 3 lines of data file:\n")
      print(head(dat,n=3L))
    }
    
   if(output==1) return(as.matrix(dat))
   if(output==2) return(data.frame(dat))

### END function readMatrix
}