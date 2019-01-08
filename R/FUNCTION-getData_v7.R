### This function is a component of astrochron: An R Package for Astrochronology
### Copyright (C) 2018 Stephen R. Meyers
###
###########################################################################
### getData: download data from astrochon server. 
###          (SRM: October 21, 2015; April 13, 2018; June 3, 2018; 
###                November 27, 2018; December 14, 2018, January 7, 2019)
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

        if(dat=="607-18O")
        {
          cat(" * Downloading benthic foraminifera d18O data series from Site 607\n\n")
          cat("   Please cite: Lisiecki, L.E., Raymo, M.E., 2005, A Pliocene-Pleistocene\n") 
          cat("   stack of 57 globally distributed d18O records, Paleoceanography, 20,\n")
          cat("   PA1003, doi:10.1029/2004PA001071.\n")
          cat("   [please also see additional references therein]\n") 
          download.file("http://www.geology.wisc.edu/~smeyers/astrochron/607-18O.txt.bz2",tempDat)
        }

        if(dat=="AEB-18O")
        {
          cat(" * Downloading planktonic foraminifera d18O data series from AEB\n\n")
          cat("   Please cite: van der Laan, E.,Hilgen, F.J., Lourens, L.J., de Kaenel,\n") 
          cat("   E.,Gaboardi, S., Iaccarino, S., 2012, Astronomical forcing of Northwest\n")
          cat("   African climate and glacial history during the late Messinian (6.5-5.5 Ma):\n")
          cat("   Palaeogeography, Palaeoclimatology, Palaeoecology, 313-314, 107-126.\n") 
          download.file("http://www.geology.wisc.edu/~smeyers/astrochron/AEB-18O.txt.bz2",tempDat)
        }

        if(dat=="Newark-rank")
        {
          cat(" * Downloading depth rank data from the Newark Basin\n\n")
          cat("   Please cite: Olsen, P.E., Kent, D.V., 1999, Long-period \n") 
          cat("   Milankovitch cycles from the Late Triassic and Early Jurassic\n")
          cat("   of Eastern North America and their implications for the calibration\n")
          cat("   of the early Mesozoic time-scale and the long-term behavior of\n")
          cat("   the planets: Phil. Trans., 357, 1761-1787.\n") 
          download.file("http://www.geology.wisc.edu/~smeyers/astrochron/Newark-rank.txt.bz2",tempDat)
        }

        if(dat=="CDL-rank")
        {
          cat(" * Downloading CDL rank data from the Latemar Limestone\n\n")
          cat("   Please cite: Preto, N., Hinnov, L.A., Hardie, L.A., De Zanche, V., 2001,\n") 
          cat("   Middle Triassic orbital signature recorded in the shallow-marine Latemar\n")
          cat("   carbonate buildup (Dolomites, Italy): Geology, 29, 1123-1126.\n")
          download.file("http://www.geology.wisc.edu/~smeyers/astrochron/CDL-rank.txt.bz2",tempDat)
        }

        if(dat=="DVCP2017-18O")
        {
          cat(" * Downloading d18O megasplice data\n\n")
          cat("   Please cite: De Vleeschouwer, David, Vahlenkamp, M., Crucifix, M.,\n") 
          cat("   Paelike, H., 2017, Alternating Southern and Northern Hemisphere climate\n")
          cat("   response to astronomical forcing during the past 35 m.y.: Geology, 45,\n")
          cat("   375-378.\n")
          download.file("http://www.geology.wisc.edu/~smeyers/astrochron/DVCP2017-18O.txt.bz2",tempDat)
        }


        cat(" * Decompressing\n")
        dat <- read.table(bzfile(tempDat),header=T)
  
        return(dat)

### END function getData
}
