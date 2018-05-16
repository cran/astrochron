### This function is a component of astrochron: An R Package for Astrochronology
### Copyright (C) 2018 Stephen R. Meyers
###
###########################################################################
### getData: download data from astrochon server. 
###          (SRM: October 21, 2015; April 13, 2018)
###
###########################################################################


getData <- function (dat="1262-a*")
{
        
        tempDat <- tempfile()
# download data
        if(dat=="926B-18O")
        {
          cat(" * Downloading adjusted benthic foraminifera oxygen isotope data from ODP Site 926B\n\n")
          cat("   Please cite: Paelike, H., Frazier, J., Zachos, J.C., 2006,\n") 
          cat("   Extended orbitally forced palaeoclimatic records from the  \n")
          cat("   equatorial Atlantic Ceara Rise: Quaternary Science Reviews, 25(23-24), 3138-3149\n")
          download.file("http://www.geology.wisc.edu/~smeyers/astrochron/926B-18O.txt.bz2",tempDat)
        }
  
        if(dat=="1262-a*")
        {
          cat(" * Downloading a* color data from ODP Site 1262\n\n")
          cat("   Please cite: Zachos, J.C., Kroon, D., Blum, P., et al., 2004,\n") 
          cat("   Proc. ODP, Init. Repts., 208: College Station, TX (Ocean Drilling\n")
          cat("   Program). doi:10.2973/odp.proc.ir.208.2004\n")
          download.file("http://www.geology.wisc.edu/~smeyers/astrochron/1262_a.txt.bz2",tempDat)
        }

        if(dat=="graptolite")
        {
          cat(" * Downloading graptolite Hidden Markov Model turnover probability series\n\n")
          cat("   Please cite: Crampton, J.S., Meyers, S.R., Cooper, R.A., et al., 2018\n") 
          cat("   Pacing of Paleozoic Macroevolutionary Rates by Milankovitch Grand\n")
          cat("   Cycles, PNAS.\n")
          download.file("http://www.geology.wisc.edu/~smeyers/astrochron/graptolite.txt.bz2",tempDat)
        }

        if(dat=="Xiamaling-CuAl")
        {
          cat(" * Downloading Cu/Al data series from the Xiamaling Formation\n\n")
          cat("   Please cite: Zhang, S., Wang, X., Hammarlund, E.U., et al., 2015\n") 
          cat("   Orbital forcing of climate 1.4 billion years ago, PNAS,\n")
          cat("   www.pnas.org/cgi/doi/10.1073/pnas.1502239112.\n")
          download.file("http://www.geology.wisc.edu/~smeyers/astrochron/Xiamaling-CuAl.txt.bz2",tempDat)
        }

        cat(" * Decompressing\n")
        dat <- read.table(bzfile(tempDat),header=T)
  
        return(dat)

### END function getData
}
