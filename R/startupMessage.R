.onAttach <- function(lib, pkg) 
 {
   if(interactive()) packageStartupMessage('Welcome to astrochron v1.5 (2025-04-28)\n',' Type ?astrochron to learn more\n',domain=NA, appendLF=TRUE)
   if(!interactive()) packageStartupMessage('Welcome to astrochron v1.5 (2025-04-28)\n',domain=NA, appendLF=TRUE)
 } 