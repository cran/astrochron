### This function is a component of astrochron: An R Package for Astrochronology
### Copyright (C) 2024 Stephen R. Meyers
###
###########################################################################
### eAsm function - (SRM: January 13, 2014, January 29, 2015; May 18, 2015;
###                       September 3, 2021; May 27, 2024)
###
### wrapper to conduct Evolutive ASM analysis using FORTRAN code
###
### NOTE: this code uses asm1.8_R.f
###########################################################################

eAsm <- function (spec,siglevel=0.9,target,fper=NULL,rayleigh,nyquist,sedmin=1,sedmax=5,numsed=50,linLog=1,iter=100000,ydir=1,palette=2,output=4,genplot=F)
{

cat("\n *****************************************************************\n")

cat("\nNOTE: The 'asm' and 'eAsm' functions are not available in astrochron version 1.3,
due to a segmentation fault error associated with Debian Linux and Fedora Linux 
clang compiliers.\n")

cat("\nIf you would like to use these functions, they work on the 10 other CRAN-supported 
compilers for Fedora Linux, Debian Linux, Mac and Windows. The functions are available in astrochron version 1.2.\n")

cat("\nastrochron 1.2 is archived on CRAN here: https://cran.r-project.org/src/contrib/Archive/astrochron/ \n")

cat("\nPlease note that 'asm' (Meyers and Sageman, 2007) was the predecessor to 'timeOpt' 
(Meyers, 2015), which provides a more comprehensive astrochronologic calibration and testing approach.\n")

cat("\n *****************************************************************")

#### END function eAsm
}
