### This function is a component of astrochron: An R Package for Astrochronology
### Copyright (C) 2022 Stephen R. Meyers
###
###########################################################################
### function calcPeriods- (SRM: June 16, 2017; August 22, 2017; 
###                             January 5, 2022; August 2, 2022)
# function to calculate periods in ka, given g, s and k (s is optional)
#
###########################################################################
  
calcPeriods <- function(g,s=NULL,k,output=1)
 {

  if(dim(data.frame(g))[2] > 1)  ipts <- dim(data.frame(g))[1]
  if(dim(data.frame(g))[2] == 1)  ipts = 1
  e1 <- double(ipts)
  e2 <- double(ipts)
  e3 <- double(ipts)
  e4 <- double(ipts)
  e5 <- double(ipts)
  
if(!is.null(s))
 {  
   o1 <- double(ipts)
   o2 <- double(ipts)
   o3 <- double(ipts)
   o4 <- double(ipts)
   o5 <- double(ipts)
   o6 <- double(ipts)
 }
  
  p1 <- double(ipts)
  p2 <- double(ipts)
  p3 <- double(ipts)
  p4 <- double(ipts)
  p5 <- double(ipts)

if(ipts > 1)
 {
# eccentricity periods in ka
# in order of relative amplitude
# 405.6 ka: g2-g5
  e1 = 1296/(g[,2]-g[,5])
# 94.9 ka: g4-g5
  e2 = 1296/(g[,4]-g[,5])
# 123.9 ka: g4-g2
  e3= 1296/(g[,4]-g[,2])
# 98.9 ka: g3-g5
  e4= 1296/(g[,3]-g[,5])
# 130.7 ka: g3-g2
  e5= 1296/(g[,3]-g[,2])

if(!is.null(s))
 {    
# obliquity periods in ka
# in order of relative amplitude  
# 41.0 ka: k+s3
  o1= 1296/(k+s[,3])
# 39.7 ka: k+s4
  o2= 1296/(k+s[,4])
# 40.3 ka: k+s3+g4-g3
  o3= 1296/(k+s[,3]+g[,4]-g[,3])
# 53.7 ka: k+s6
  o4= 1296/(k+s[,6]) 
# 41.7 ka: k+s3-g4+g3
  o5= 1296/(k+s[,3]-g[,4]+g[,3])
# 28.9 ka: k+s1
  o6= 1296/(k+s[,1])  
 }   
  
# precession periods in ka
# in order of relative amplitude
# 23.7 ka: k+g5
  p1= 1296/(k+g[,5])
# 22.4 ka: k+g2
  p2= 1296/(k+g[,2])
# 18.9 ka: k+g4
  p3= 1296/(k+g[,4])
# 19.1 ka: k+g3
  p4= 1296/(k+g[,3]) 
# 23.1 ka: k+g1
  p5= 1296/(k+g[,1])
 } 
  
if(ipts == 1)
 {
# eccentricity periods in ka
# in order of relative amplitude
# 405.6 ka: g2-g5
  e1 = 1296/(g[2]-g[5])
# 94.9 ka: g4-g5
  e2 = 1296/(g[4]-g[5])
# 123.9 ka: g4-g2
  e3= 1296/(g[4]-g[2])
# 98.9 ka: g3-g5
  e4= 1296/(g[3]-g[5])
# 130.7 ka: g3-g2
  e5= 1296/(g[3]-g[2])
  
if(!is.null(s))
 {    
# obliquity periods in ka
# in order of relative amplitude  
# 41.0 ka: k+s3
  o1= 1296/(k+s[3])
# 39.7 ka: k+s4
  o2= 1296/(k+s[4])
# 40.3 ka: k+s3+g4-g3
  o3= 1296/(k+s[3]+g[4]-g[3])
# 53.7 ka: k+s6
  o4= 1296/(k+s[6]) 
# 41.7 ka: k+s3-g4+g3
  o5= 1296/(k+s[3]-g[4]+g[3])
# 28.9 ka: k+s1
  o6= 1296/(k+s[1])  
 }     

# precession periods in ka
# in order of relative amplitude
# 23.7 ka: k+g5
  p1= 1296/(k+g[5])
# 22.4 ka: k+g2
  p2= 1296/(k+g[2])
# 18.9 ka: k+g4
  p3= 1296/(k+g[4])
# 19.1 ka: k+g3
  p4= 1296/(k+g[3]) 
# 23.1 ka: k+g1
  p5= 1296/(k+g[1])
 }  
  
# output from longest to shortest period  
  if(is.null(s)) 
   {
     periodsOut=cbind(e1,e5,e3,e4,e2,p1,p5,p2,p4,p3)
     colnames(periodsOut) <- c("e1","e5","e3","e4","e2","p1","p5","p2","p4","p3")
   } 

# output from longest to shortest period  
  if(!is.null(s)) 
   {
     periodsOut=cbind(e1,e5,e3,e4,e2,o4,o5,o1,o3,o2,o6,p1,p5,p2,p4,p3)
     colnames(periodsOut) <- c("e1","e5","e3","e4","e2","o4","o5","o1","o3","o2","o6","p1","p5","p2","p4","p3")
   }   
    
  if(output==1) return(data.frame(periodsOut)) 
  if(output==2) return(as.numeric(periodsOut)) 
 } 