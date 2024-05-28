.onAttach <- function(lib, pkg) 
 {
   if(interactive()) packageStartupMessage('Welcome to astrochron v1.3 (2024-05-27)\n',' Type ?astrochron to learn more\n',domain=NA, appendLF=TRUE)
   if(!interactive()) packageStartupMessage('Welcome to astrochron v1.3 (2024-05-27)\n',domain=NA, appendLF=TRUE)
 } 