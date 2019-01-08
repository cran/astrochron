.onAttach <- function(lib, pkg) 
 {
   if(interactive()) packageStartupMessage('Welcome to astrochron v0.9 (2019-01-08)\n',' Type ?astrochron to learn more\n',domain=NA, appendLF=TRUE)
   if(!interactive()) packageStartupMessage('Welcome to astrochron v0.9 (2019-01-08)\n',domain=NA, appendLF=TRUE)
 } 