.onAttach <- function(lib, pkg) 
 {
   if(interactive()) packageStartupMessage('Welcome to astrochron v1.1 (2021-11-05)\n',' Type ?astrochron to learn more\n',domain=NA, appendLF=TRUE)
   if(!interactive()) packageStartupMessage('Welcome to astrochron v1.0 (2021-11-05)\n',domain=NA, appendLF=TRUE)
 } 