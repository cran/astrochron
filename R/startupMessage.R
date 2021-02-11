.onAttach <- function(lib, pkg) 
 {
   if(interactive()) packageStartupMessage('Welcome to astrochron v1.0 (2021-02-10)\n',' Type ?astrochron to learn more\n',domain=NA, appendLF=TRUE)
   if(!interactive()) packageStartupMessage('Welcome to astrochron v1.0 (2021-02-10)\n',domain=NA, appendLF=TRUE)
 } 