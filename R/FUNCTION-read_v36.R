### This function is a component of astrochron: An R Package for Astrochronology
### Copyright (C) 2023 Stephen R. Meyers
###
###########################################################################
### Read function - (SRM: January 24, 2012; updated January 31, 2012; 
###                  March 5, 2012; April 20, 2012; June 30, 2012; 
###                  December 2, 2012; March 11, 2013; March 12, 2013; May 1-2, 2013;
###                  May 17, 2013; June 5, 2013; June 28, 2013; July 26-27, 2013;
###                  Nov. 27, 2013; January 14, 2014; June 25, 2014; January 22, 2015;
###                  August 6, 2015; September 11, 2015; September 16, 2015; 
###                  October 8, 2015; November 19, 2016; March 20, 2017; July 26, 2017;
###                  December 15, 2017; July 4, 2018; August 13, 2018; January 14, 2021;
###                  June 22, 2021; November 7, 2022; July 27, 2023)
###
### Read a time series file. 
###########################################################################

read <- function (file=NULL,d=1,h="auto",skip=0,srt=T,ave=T,check=T,genplot=T,verbose=T)

{
 if(verbose)
  {
   cat("\n----- READ STRATIGRAPHIC SERIES FROM DATA FILE -----\n")
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
       xlab="Location"
       ylab="Value"
     }

  if(h=="yes")
     {
       if(d == 0) {dat <- read.table (filen,header=T,skip=skip)}
       if(d == 1) {dat <- read.table (filen,header=T,sep=",",skip=skip)}
       if(d == 2) {dat <- read.table (filen,header=T,sep=";",skip=skip)}
       xlab=names(dat[1])
       ylab=names(dat[2])
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
           xlab="Location"
           ylab="Value"
         }
# column titles        
        if (titles )  
         {
           if(verbose) cat("\n * Column titles/headers detected\n")
           dat <- read.table (filen,header=T,skip=skip)
           xlab=names(dat[1])
           ylab=names(dat[2]) 
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
           ylab="Value"         
         }
# column titles        
        if (titles)  
         {
          if(verbose) cat("\n * Column titles/headers detected\n")
          dat <- read.table (filen,header=T,sep=",",skip=skip)
          xlab=names(dat[1])
          ylab=names(dat[2])    
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
           ylab="Value"         
         }
# column titles        
        if (titles)  
         {
          if(verbose) cat("\n * Column titles/headers detected\n")
          dat <- read.table (filen,header=T,sep=";",skip=skip)
          xlab=names(dat[1])
          ylab=names(dat[2])    
         }
      }

   npts <- length(dat[,1]) 
   if(verbose) cat("\n * Number of stratigraphic samples (rows)=", npts,"\n")

   cols=length(dat[1,])
   if(verbose) cat(" * Number of variables (columns)=",cols-1," (excluding depth/height/time)\n")

# note: if check is disabled sorting and averaging are also disabled 
   if(!check && verbose) cat("\n  WARNING: data checking disabled. Sorting, duplicate averaging and empty entry removal are also disabled.\n")
   if(check)
    {
    
### if we have more than two columns: 
   if (cols > 2)
    {
# check to see if any of the rows are all NA entries
      delRow <- logical(npts)
      for (i in 1:npts) delRow[i]=all(is.na(dat[i,]))
# delete rows that have all NA entries  
      if(any(delRow))
        { 
          dat <- dat[!delRow,]
          npts <- length(dat[,1]) 
          if(verbose) cat("\n * Some rows contain all NA entries, and will be removed\n")
          if(verbose) cat(" * New number of rows=", npts,"\n")
        }
# check to see if any of the columns are all NA entries
      delCol<-logical(cols)
      for (i in 1:cols) delCol[i]=all(is.na(dat[,i]))
# delete rows that have all NA entries  
      if(any(delCol))
        { 
          dat<-dat[,!delCol]
          cols <- length(dat[1,]) 
          if(verbose) cat("\n * Some columns contain all NA entries, and will be removed\n")
          if(verbose) cat(" * New number of columns=", cols-1," (excluding depth/height/time)\n")
        }
    }

### if we have two columns (time, variable), then sort, remove if NA in either column, and 
###  average duplicates if requested  
   if (cols == 2 && srt)
    {
       if(verbose) cat("\n * Sorting data into increasing depth/height/time order.\n")
       if(verbose) cat("   Will remove empty entries (from either column).\n")
### Remove NAS entries in second column (includes those listed as 'NA')
       dat=subset(dat,!is.na(dat[,2]))
### sort to ensure increasing depth/height/time. Will remove NAS values (in depth column)
       dat <- dat[order(dat[,1],na.last=NA,decreasing=F),]
       npts <- length(dat[,1])
       if(verbose) cat(" * Number of samples (rows) post-sorting=", npts,"\n")

### function dup: average duplicates/triplicates/etc.
dup <- function (ipts,x,y)
 {
    F_dat = .Fortran('dupmean_r',
    
    ipts=as.integer(ipts),x=as.double(x),
    y=as.double(y),npts=integer(1),xout=double(ipts),yout=double(ipts))

# return the results
    return(cbind(F_dat$xout[1:F_dat$npts],F_dat$yout[1:F_dat$npts]))
  }
### END function dup

### if there are duplicate values, average them
      t1<-dat[1:(npts-1),1]
      t2<-dat[2:(npts),1]
      dt=t2-t1
      mindt=min(dt) 
# save a copy of dat before averaging   
      dat1=dat    
       if(mindt < 1e-9)
        {
         if(verbose) cat(" * Duplicates found\n")
         if( ave && ( inherits(dat[,2],"numeric") || inherits(dat[,2],"integer") ) )
           {  
             if(verbose) cat(" * Duplicates values will be averaged.\n")
### call to Fortran routine for quick duplicate averaging
             dat2 <- dup(npts,dat[,1],dat[,2])
             dat2 <- data.frame(dat2)
             npts <- length(dat2[,1]) 
             if(verbose) cat(" * New number of samples (rows)=",npts,"\n")
             colnames(dat2)[1] <- colnames(dat[1])
             colnames(dat2)[2] <- colnames(dat[2])
             dat <- dat2
            }
         else if(verbose) cat(" * Duplicates found, but will not be averaged.\n")
        } 
    }

### if we have more than two columns, then sort and remove NA's in first column
###   note that duplicates are not averaged.
   if (cols > 2 && srt)
    {
      if(verbose)
       {
         cat("\n * Sorting data into increasing depth/height/time order.\n")
         cat("   Will remove empty entries from depth/height/time column only\n")
         cat("   (empty entries may remain in other columns).\n")
       }  
### sort to ensure increasing depth/height/time. Will remove NAS values (in depth column)
      dat <- dat[order(dat[,1],na.last=NA,decreasing=F),]
      npts <- length(dat[,1])
      if(verbose) cat("\n * Number of rows (samples) post-sorting=", npts,"\n")
### determine if there are duplicate values
      t1<-dat[1:(npts-1),1]
      t2<-dat[2:(npts),1]
      dt=t2-t1
      mindt=min(dt) 
      if(mindt < 1e-9)
       {
         if(verbose)
          {
            cat(" * Duplicates found\n")
            if(ave) cat("\n  WARNING: Cannot average duplicate values\n")
          }  
       }  
    }

# determine if any missing entries (NA) still exist (this may be the case when cols>2)
     numNA=sum(is.na(dat))
     if(numNA > 0) 
       {
         if(verbose) cat("\n  WARNING:", numNA,"empty entries are still present in your data frame.\n")
       } 

# end check section
    }
  
  
### now evaluate sampling statistics
     t1<-dat[1:(npts-1),1]
     t2<-dat[2:(npts),1]
     dt=t2-t1
     dtMin=min(dt) 
     dtMax=max(dt)
     dtMean=mean(dt)     
     dtMedian=median(dt)

     if(verbose) 
      {
        cat("\n * Mean sampling interval=", dtMean,"\n")
        cat(" * Median sampling interval=",dtMedian,"\n")
        cat(" * Maximum sampling interval=",dtMax,"\n")
        cat(" * Minimum sampling interval=", dtMin,"\n")
      }
      
if(genplot && cols==2)
 {
   if( inherits(dat[,2],"numeric") || inherits(dat[,2],"integer") )
    {
### plot data series. Note, cex is the factor by which to increase or decrease default symbol size
      par(mfrow=c(2,2))
      if( ave && ( inherits(dat[,2],"numeric") || inherits(dat[,2],"integer") ) )
       {
        plot(dat,type="l",col="gray",ylim=c(min(dat1[,2]),max(dat1[,2])),xlab=xlab,ylab=ylab,main="Stratigraphic Series")
        mtext("(red=original data, black=post-averaging)",cex=0.8)
        points(dat1,cex=0.4,col="red")
        points(dat,cex=0.3)         
       }
      else{   
         plot(dat,type="l",col="gray",xlab=xlab,ylab=ylab,main="Stratigraphic Series")
         points(dat,cex=0.3)
       }
### plot the density and the histogram together
      hist(dat[,2],freq=F,xlab=ylab,main=paste("Distribution of",ylab)); grid(); lines(density(dat[,2], bw="nrd0",na.rm=T),col="red",lwd=2)
### boxplot
      boxplot(dat[,2],ylab=ylab,main=paste("Boxplot of",ylab))
### Normal probabilty plot (Normal Q-Q Plot)
      qqnorm(dat[,2],main=paste("Normal Q-Q plot of",ylab)); qqline(dat[,2], col="red",lwd=2); grid()
    }
   if( !inherits(dat[,2],"numeric") && !inherits(dat[,2],"integer") )
    {
      if(verbose) cat("\n  WARNING: data contains non-numeric values.  Will NOT generate plots.\n")
    }
  } 
 
if(genplot && cols>2) 
  {
    nrow = ceiling(sqrt(cols-1))
    ncol = nrow
    par(mfrow = c(nrow, ncol))
    for (i in 2:cols) 
      {
        xlab = names(dat)[i]
        plotdat = subset(dat[,i], !(dat[,i] == "NA"))
        if( inherits(dat[,i],"numeric") || inherits(dat[,i],"integer") )
          {
            hist(plotdat,freq=F,xlab=xlab,main=""); lines(density(plotdat, bw="nrd0",na.rm=T),col="red",lwd=2)
          } 
        if( !inherits(dat[,i],"numeric") && !inherits(dat[,i],"integer") )
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
   
   return(data.frame(dat))

### END function read
}