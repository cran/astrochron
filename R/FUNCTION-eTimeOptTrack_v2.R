### This function is a component of astrochron: An R Package for Astrochronology
### Copyright (C) 2018 Stephen R. Meyers
###
###########################################################################
### eTimeOptTrack function - (SRM: October 5, 2016; December 7, 2017)
###
### track eTimeOpt maxmia
###########################################################################

eTimeOptTrack <- function (res,threshold=0,ydir=-1,genplot=T,verbose=T)
{

 if(verbose) cat("\n----- IDENTIFYING MAXIMUM r2 value IN EACH WINDOW -----\n")

# ensure we have a data frame
  res=data.frame(res)

# assign sedrates from first column of res
  sedrates=res[,1]

  cols=length(res)
# assign locations for each spectrum (column headers)
  loc=suppressWarnings(as.numeric(substr(colnames(res[2:cols]),start=2,stop=100)))
# for negative depth/height/time values, "-" has been changed to "."
# this will create NAs. implement modification of fix recommended by Mathieu Martinez
  neg=grepl(".",substr(colnames(res[2:cols]), start=2,stop=2),fixed=T)
  fixloc=which(neg)
  if(any(neg)) {loc[fixloc]=-1*as.numeric(substr(colnames(res[(fixloc+1)]),start=3,stop=100))}

# assign r2 results
  r2=as.matrix( res[2:length(res)] )

  numrec=length(loc)
  numsed=length(sedrates)
  if(verbose) cat("\n * Number of windows to analyze =",numrec,"\n")
  if(verbose) cat(" * Number of sedimentation rates =",numsed,"\n")

# loop over all windows 
# note: these are dimensioned to be much larger than needed
   r2Max<-double(numsed*numrec)
   sedMax<-double(numsed*numrec)
   locHold<-double(numsed*numrec)
   j=1
  for (i in 1:numrec)
    {
#     if(verbose) cat(" * PROCESSING WINDOW=",i,"; Location=",loc[i],"\n")
# find max r2
     amax <- max(r2[,i])
# identify all sedimentation rates with this r2 value. This allows for the case
#  when there are multiple equivalent maxima (unlike 'max' or 'which.max')
     imax=which(r2[,i] == amax,arr.ind=TRUE)
     if(length(imax) > 1 && verbose) cat( "**** WARNING IN WINDOW",i,"(",loc[i],"):", length(imax)," sedimentation rates have global maximum value\n")
     for (k in 1:length(imax))
       {
         r2Max[j]=r2[imax[k],i]
         sedMax[j]=sedrates[imax[k]]
         locHold[j]=loc[i]
         j=j+1
       }
# end numrec loop
    }

   out <- data.frame(cbind(locHold[1:(j-1)],sedMax[1:(j-1)],r2Max[1:(j-1)]))
   colnames(out)<-c("Location","Sedrate","r2")
   
   out <- subset(out,(out[,3] >= threshold))

  if(genplot)
   {
     par(mfrow=c(1,2))
     if(ydir == 1) ylim=c(min(out[,1]),max(out[,1]))
     if(ydir == -1) ylim=c(max(out[,1]),min(out[,1]))
     plot(out[,2],out[,1],type="b",cex=0.75,ylim=ylim,xlab="Sedimentation Rate",ylab="Location")
     plot(out[,3],out[,1],type="b",cex=0.75,ylim=ylim,xlab="r2",ylab="Location")
   }
 
   return(out)

#### END function eTimeOptTrack
}
