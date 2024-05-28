### This function is a component of astrochron: An R Package for Astrochronology
### Copyright (C) 2024 Stephen R. Meyers
###
###########################################################################
### stratPCA: perform basic principal component analysis 
###            (SRM Sept. 19-21, 2023; Oct. 5-7, 2023; May 27, 2024)
###
###########################################################################

stratPCA <- function (dat,id=TRUE,rot=0,nPC=NULL,output=0,symSize=2,genplot=1)
{

### id: first column is an ID column? T or F
### rot: choose rotation, 0= none, 1= varimax, 2= promax
### nPC: number of principal components to extract (default = # variables)
### genplot: generate plots? (1) standard plots, (2) additional plots
### output: 0 = nothing, 1 = loadings; 2 = scores

    ipts <- length(dat[,1]) 
    cat("\n * Number of stratigraphic data/observation levels=", ipts, "\n")

### remove depth column for further statistical analayses if need be
    if(id) { work <- dat[,2:ncol(dat)] }
### otherwise
    if(!id) { work <- dat }

    if(is.null(nPC)) nPC=ncol(work)
    nvar=ncol(work)
    
    if(rot == 0) { rotation="none" }
    if(rot == 1) { rotation="varimax" }
    if(rot == 2) { rotation="promax" }

# (1) caculate the correlation matrix
    cat(" * Calculating correlation matrix \n")
    correl <- cor(work)
    
  if(genplot >= 2)
   {
    dev.new(title="Correlation Plot",height=4.5,width=5)
    par(mar=c(4,5,4,3)) 
    par(fig=c(0,0.90,0,1))
    plot(-1,-1,ylim=c(1,nvar),xlim=c(1,nvar),xlab="",ylab="",xaxt="n",yaxt="n")
    pal <- colorRamp(c("blue", "white", "red"))
    for(i in 1:nvar)
     { 
      col <- rgb(pal((correl[,i] - (-1)) / 2), maxColorValue=255)
      points(rep(i,nvar),1:nvar,col=col,pch=15,cex=symSize)
     }
    axis(3, at=1:nvar, labels=colnames(correl))
    axis(2, at=1:nvar, labels=rownames(correl),las=1) 
    par(new = TRUE)
    par(mar=c(0.1,0.1,0.1,1)) 
    par(fig=c(0.85,0.95,0.2,0.8))
    image(1, seq(0,1,.01), t(seq_along(seq(0,1,.01))), col=rgb(pal(seq(0,1,.01)),maxColorValue=255), axes=FALSE,xlab="",ylab="")
    text(1,0.95,"1",cex=1.1,col="white")
    text(1,0.05,"-1",cex=1.1,col="white")
   }  

# (2) using eigen instead of svd, as testing indicates that it is consistent with SAS approach. also consistent with S-PLUS.
    cat(" * Performing eigenvalue-eigenvector value decomposition \n\n")
    eigens <- eigen(correl)

# (3) extract eigenvalues
    pca.eig <- eigens$values

# (4) calculate loadings using matrix multiplication
###     eigenvectors * sqrt(eigenvalues)
    pca.load <- eigens$vectors %*% sqrt(diag(pca.eig))
    rownames(pca.load) <- rownames(correl)  
   
### calculate scores: multiply data matrix by loading
###   matrix (matrix multiplication), as in Davis (p.504) and Iacobucci notes (Northwestern Univ.).
    score = as.matrix(work) %*% pca.load
   
#### output unrotated results   
    cat("\n INITIAL COMMUNALITY ESTIMATES: \n") 
    cat(rep(1,ncol(work)))
    cat("\n\n EIGENVALUES of CORRELATION MATRIX: \n") 
    print(pca.eig)

### show scree plot for unrotated results
if(genplot >= 1)
 {
    dev.new(title="Eigenvalues: scree plot (No Rotation)",height=7.5,width=4)
    par(mfrow=c(2,1))
    par(mar=c(4,3,2,3))
    plot(pca.eig,xlab="",ylab="",type="b")
    mtext("Eigenvalue",side=2,line=1.9,cex=0.9)
    mtext("Principal Component",side=1,line=1.9,cex=0.9)
    abline(h=1,lty=2)
    par(new=T)
    plot(100*pca.eig/length(work),xlab="",ylab="",type="b",axes=F,lwd=2)
    axis(4, ylim=c(0,max(100*pca.eig/length(work)),lwd=1,col="black"))
    mtext("% Variance",side=4,line=1.9,cex=0.9) 
    plot(log(pca.eig),xlab="",ylab="",type="b",lwd=2)
    mtext("Log(Eigenvalue)",side=2,line=1.9,cex=0.9)
    mtext("Principal Component",side=1,line=1.9,cex=0.9)
 }

if(genplot >= 1 && rot == 0)
 {
# plot loadings 
   dev.new(title=paste("Loadings: Number of PCs=",nPC,", Rotation=None"),height=4.5,width=5)
   par(mar=c(4,5,4,3)) 
   par(fig=c(0,0.90,0,1))
   plot(-1,-1,ylim=c(1,nvar),xlim=c(1,nPC),xlab="",ylab="",xaxt="n",yaxt="n")
   mtext("Principal Component",line=2)
   pal <- colorRamp(c("blue", "white", "red"))
   for(i in 1:nPC)
    { 
     col <- rgb(pal((pca.load[,i] - (-1)) / 2), maxColorValue=255)
     points(rep(i,nvar),1:nvar,col=col,pch=15,cex=symSize)
    } 
    axis(3, at=1:nPC, labels=colnames(pca.load))
    axis(2, at=1:nvar, labels=rownames(pca.load),las=1)
    par(new = TRUE)
    par(mar=c(0.1,0.1,0.1,1)) 
    par(fig=c(0.85,0.95,0.2,0.8))
    image(1, seq(0,1,.01), t(seq_along(seq(0,1,.01))), col=rgb(pal(seq(0,1,.01)),maxColorValue=255), axes=FALSE,xlab="",ylab="")
    text(1,0.95,"1",cex=1.1,col="white")
    text(1,0.05,"-1",cex=1.1,col="white")

### show multiple score plots if desired (up to 10), on new device
    dev.new(title=paste("Scores vs. Depth: Number of PCs=",nPC,", Rotation=None"))
    if(nPC<=5) {nrows=1;ncols=nPC}
    if(nPC>5 ) {nrows=2;ncols=5}
    if(nPC<=10) maxPC=nPC
    if(nPC>10) maxPC=10
    par(mfrow=c(nrows,ncols))
    par(mar=c(4,3,2,1))
### loop over PCs of interest
    for (i in 1:maxPC)
     {
       plot(score[,i],dat[,1]*-1,type="l", ylab="", xlab="",cex=0.5,bty="n",lwd=0.5,main=paste("PC",i))
       mtext("PC score",side=1,line=2, cex=0.7)
     }
 }

### show loading and scores scatter plots if desired
if(genplot >= 2 && nPC == ncol(work))
 {
### crossplots of PC loadings and scores (up to 5 components)
    if(nPC==2) {nrows=1; ncols=1; setht=8; setwd=8}
    if(nPC==3) {nrows=2; ncols=2; setht=8; setwd=8}
    if(nPC==4) {nrows=2; ncols=3; setht=7; setwd=10}
    if(nPC>=5) {nrows=3; ncols=4; setht=7; setwd=10}
    if(nPC<=5) maxPC=nPC
    if(nPC>5) maxPC=5
    dev.new(title=paste("Loadings vs. Loadings: Number of PCs=",nPC,", No Rotation"),height=setht,width=setwd)
    par(mfrow=c(nrows,ncols))
    par(mar=c(4,4,2,2))
### loop over PCs of interest
    for (i in 1:(maxPC-1))
    {
     for (ii in (i+1):(maxPC))
      {
        plot(pca.load[,i],pca.load[,ii],xlim=c(-1,1),ylim=c(-1,1),xlab="",ylab="",bty="n",font.lab=2)
        mtext(paste("PC",i,"loading"),side=1,line=2.5,cex=0.7,font=2)
        mtext(paste("PC",ii,"loading"),side=2,line=2.5,cex=0.7,font=2)
        text(pca.load[,i],pca.load[,ii], label = rownames(pca.load), col="red",cex=0.75)
      }
     } 
    dev.new(title=paste("Scores vs. Scores: Number of PCs=",nPC,", No Rotation"),height=setht,width=setwd)
    par(mfrow=c(nrows,ncols))
    par(mar=c(4,4,2,2))
    for (i in 1:(maxPC-1))
    {
     for (ii in (i+1):(maxPC))
      {
        plot(score[,i],score[,ii],xlab="",ylab="",bty="n",font.lab=2)
        mtext(paste("PC",i,"score"),side=1,line=2.5,cex=0.7,font=2)
        mtext(paste("PC",ii,"score"),side=2,line=2.5,cex=0.7,font=2)
      }
     } 
  }
### summary of initial PCA

    cat("\n PRE-ROTATION PC LOADINGS: \n")
    print(round(pca.load[,1:nPC],digits=3),cutoff=0)
 
    cat("\n ->",nPC,"PCs account for", 100*sum(pca.eig[1:nPC])/length(work),"% of the variance\n")

if(rot > 0 && nPC==ncol(work)) 
 {
   cat("\n **** WARNING: Rotation will not be performed unless nPC < number of input variables." )
   cat("\n               Please reduce nPC, and rerun. \n" )
   rot=0
 }
   

####################################################
# allow rotation if extracting nPC < # variables
#################################################### 
   
if(nPC<(ncol(work)) && rot > 0)
{

    cat("\n * Performing rotation on loadings \n")
### extract nPC loadings
     pca.sub <- pca.load[,1:nPC]
### conduct a VARIMAX ROTATION- NOTE: eps set here to match SAS results
     if(rot==1)  pca.rot.load <- varimax(pca.sub,eps=1e-12)$loadings
### conduct a PROMAX ROTATION- NOTE: power used for target (m) should be 2-4  
     if(rot==2) pca.rot.load <- promax(pca.sub,m=4)$loadings
   
### summary of final FA
    cat("\n POST-ROTATION COMMUNALITY ESTIMATES: \n")
    print(rowSums(pca.rot.load^2))

### PROVIDE A MORE COMPLETE LIST OF LOADINGS FOR OUTPUT
   comp.load <- matrix(nrow=ncol(work),ncol=nPC)
   for (i in 1:nPC)
    {
     comp.load[,i]  <- pca.rot.load[,i]
    }
   colnames(comp.load) <- colnames(pca.rot.load)
   rownames(comp.load) <- rownames(pca.rot.load)

   cat("\n POST-ROTATION PC LOADINGS:")
   print(round(pca.rot.load, digits=3),cutoff=0)
   
# caculate new eigenvalues for rotated matrix (sum of the squares of PC loadings)
   rot.eig <- double(nPC)
   for (i in 1:nPC)
    {
     rot.eig[i] <- sum(pca.rot.load[,i]^2)
    }
   

### calculate scores: multiply data matrix by loading
###   matrix (matrix multiplication), as in Davis (p.504) and Iacobucci notes (Northwestern Univ.).
    score.rot = as.matrix(work) %*% pca.rot.load  # Dec. 12, 2019, keep original order

  if(rot==1) { minVal=-1 ; maxVal=1; rangeVal=2}
#  note that oblique rotations can yield loadings >1 or <-1
#  center on the color scale even when values are >1 or <-1
   if(rot==2) { minVal= min (min(-1,min(comp.load)), -1*max(1,max(comp.load))) ; maxVal= max(max(comp.load),-1*minVal,1) ; rangeVal=maxVal-minVal }
   
if(genplot >= 1)
 {
# plot loadings 
   dev.new(title=paste("Loadings: Number of PCs=",nPC,", Rotation=",rotation),height=4.5,width=5)
   par(mar=c(4,5,4,3)) 
   par(fig=c(0,0.90,0,1))
   plot(-1,-1,ylim=c(1,nvar),xlim=c(1,nPC),xlab="",ylab="",xaxt="n",yaxt="n")
   mtext("Principal Component",line=2)
   pal <- colorRamp(c("blue", "white", "red"))
  for(i in 1:nPC)
   { 
    col <- rgb(pal((comp.load[,i] - (minVal)) / rangeVal), maxColorValue=255)
    points(rep(i,nvar),1:nvar,col=col,pch=15,cex=symSize)
   } 
   axis(3, at=1:nPC, labels=colnames(comp.load))
   axis(2, at=1:nvar, labels=rownames(comp.load),las=1)
   par(new = TRUE)
   par(mar=c(0.1,0.1,0.1,1)) 
   par(fig=c(0.85,0.95,0.2,0.8))
   image(1, seq(0,1,.01), t(seq_along(seq(0,1,.01))), col=rgb(pal(seq(0,1,.01)),maxColorValue=255), axes=FALSE,xlab="",ylab="")
   text(1,0.95,round(maxVal,1),cex=1.1,col="white")
   text(1,0.05,round(minVal,1),cex=1.1,col="white")
  
### show multiple score plots if desired (up to 10), on new device  
    dev.new(title=paste("Scores vs. Depth: Number of PCs=",nPC,", Rotation=",rotation))
    if(nPC<=5) {nrows=1;ncols=nPC}
    if(nPC>5 ) {nrows=2;ncols=5}
    if(nPC<=10) maxPC=nPC
    if(nPC>10) maxPC=10
    par(mfrow=c(nrows,ncols))
    par(mar=c(4,3,2,1))
### loop over PCs of interest
    for (i in 1:maxPC)
     {
       plot(score.rot[,i],dat[,1]*-1,type="l", ylab="", xlab="",cex=0.5,bty="n",lwd=0.5,main=paste("PC",i))
       mtext("PC score",side=1,line=2, cex=0.7)
     }
 }

### show loading scatter plots if desired
if(genplot >= 2)
 {
### crossplots of PC loadings (up to 5 components)
    if(nPC==2) {nrows=1; ncols=1; setht=8; setwd=8}
    if(nPC==3) {nrows=2; ncols=2; setht=8; setwd=8}
    if(nPC==4) {nrows=2; ncols=3; setht=7; setwd=10}
    if(nPC>=5) {nrows=3; ncols=4; setht=7; setwd=10}
    if(nPC<=5) maxPC=nPC
    if(nPC>5) maxPC=5
    dev.new(title=paste("Loadings vs. Loadings: Number of PCs=",nPC,", Rotation=",rotation),height=setht,width=setwd)
    par(mfrow=c(nrows,ncols))
    par(mar=c(4,4,2,2))
### loop over PCs of interest
    for (i in 1:(maxPC-1))
     {
     for (ii in (i+1):(maxPC))
       {     
         plot(pca.rot.load[,i],pca.rot.load[,ii],xlim=c(-1,1),ylim=c(-1,1),xlab="",ylab="",font.lab=2)
         mtext(paste("PC",i,"loading"),side=1,line=2.5,cex=0.7,font=2)
         mtext(paste("PC",ii,"loading"),side=2,line=2.5,cex=0.7,font=2)
         text(pca.rot.load[,i],pca.rot.load[,ii], label = rownames(pca.rot.load), col="red",cex=0.5)
       }
     } 
    dev.new(title=paste("Scores vs. Scores: Number of PCs=",nPC,", Rotation=",rotation),height=setht,width=setwd)
    par(mfrow=c(nrows,ncols))
    par(mar=c(4,4,2,2))
    for (i in 1:(maxPC-1))
    {
     for (ii in (i+1):(maxPC))
      {
        plot(score.rot[,i],score.rot[,ii],xlab="",ylab="",bty="n",font.lab=2)
        mtext(paste("PC",i,"score"),side=1,line=2.5,cex=0.7,font=2)
        mtext(paste("PC",ii,"score"),side=2,line=2.5,cex=0.7,font=2)
      }
     }

  }   
# end if(nPC<(ncol(work)))
}

if(output==1 && rot >0) return(comp.load)
if(output==1 && rot ==0) return(data.frame(pca.load))

if(output==2 && rot >0) return(data.frame(cbind(dat[,1],score.rot)))
if(output==2 && rot ==0) return(data.frame(cbind(dat[,1],score)))

  
### END function stratPCA
}
