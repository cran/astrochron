### This function is a component of astrochron: An R Package for Astrochronology
### Copyright (C) 2023 Stephen R. Meyers
###
###########################################################################
### function insoDiff : calculate difference between two insolation maps using 
###                    Laskar et al. (2004) solution, with palinsol 
###                    (SRM: November 7-9, 2022; July 22, 2023)
###
###########################################################################


insoDiff <- function (t1,t2,S0=1365,verbose=T) 
{
  if(verbose) cat("\n----- CALCULATE INSOLATION MAP USING LASKAR ET AL. (2004) SOLUTION -----\n")

# error checking
# restrict to 21 million years into future (>= -21000). 
  if(t1 < (-21000))
   {
    if(verbose) cat("\n**** WARNING: This solution can only be calculated to 21000 kiloyears in the future. Resetting t1 to -20000.\n")
    t1=-21000
   }   
# restrict to 51 million years into past (<= 51000).
  if(t1 > 51000)
   {
    if(verbose) cat("\n**** WARNING: This solution can only be calculated to 51000 kiloyears in the past. Resetting t1 to 50000.\n")
    t1=51000
   }    
# restrict to 21 million years into future (>= -21000). 
  if(t2 < (-21000))
   {
    if(verbose) cat("\n**** WARNING: This solution can only be calculated to 21000 kiloyears in the future. Resetting t2 to -20000.\n")
    t2=-21000
   }   
# restrict to 51 million years into past (<= 51000).
  if(t2 > 51000)
   {
    if(verbose) cat("\n**** WARNING: This solution can only be calculated to 51000 kiloyears in the past. Resetting t2 to 50000.\n")
    t2=51000
   }  

# convert time from kyr to years 
   M1 <- Milankovitch(la04(t1*1000))
   M2 <- Milankovitch(la04(t2*1000)) 
   Mdiff=M1-M2
   dev.new(height=3.5,width=9)
   par(mfrow = c(1, 3))
   plot(M1,plot=contour,main=paste(t1,"ka"))
   plot(M2,plot=contour,main=paste(t2,"ka"))     
   plot(Mdiff,plot=contour,main=paste(t1,"-",t2,"ka"))
   
### END function insoDiff
}

