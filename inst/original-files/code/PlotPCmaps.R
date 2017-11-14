library(sp)
library(rgdal)
library(RColorBrewer)

source("prepSpatial.R")

# Load PC data
#base="GSKd3hm.Swiss.LD_0.8.exLD-PCA"
#base="mergedSub.EuropeW.LD_0.8.exLD-PCA"
#base="mergedSub.EuroThinFinal.LD_0.8.exLD.out0-PCA"
#base="GSK_08_24_01.EuroThinFinal.LD_0.8.exLD.out0-PCA"
#base="ZoltanEvenSample-PCA"
base="POPRES_08_24_01.EuroThinFinal.LD_0.8.exLD.out0-PCA"


#filename=paste("../pca_analysis/",base,"/",base,".eigs",sep="")
filename=paste(base,".eigs",sep="")

PCA=read.table(filename)
#PCA=read.table("../pca_analysis/mergedSub.EuroNeat.LD_0.8.exLD.out0-PCA/mergedSub.EuroNeat.LD_0.8.exLD.out0-PCA.eigs")
#PCA=read.table("../pca_analysis/mergedSub.EuroNeat.topPC1PC2.LD_0.8ex.LD-PCA/mergedSub.EuroNeat.topPC1PC2.LD_0.8ex.LD-PCA.eigs")
names(PCA)=c("FamID","ID",paste("PC",1:(dim(PCA)[2]-2),sep=""))
#PCA=PCA[!is.element(demodata[as.character(PCA$ID),"Wave"],c(6,7)),]

cex.text=1

matchord=match(PCA$ID,ddata$POPRES_ID)
rlabels=ddata$tagGPnodiscord[matchord]
clabels=collabelmap[matchord]
plabels=poslabel[matchord]
alabels=abbrvlabel[matchord]
grplabels=grplabels[matchord]
maplabel=maplabel[matchord]

# See demo.biplot.R in this repository for a working example of this
# function when map = TRUE, rotate = TRUE and useCols = TRUE. -Peter
plotPCbi=function(x=1,y=2,sx=1,sy=1,map=TRUE,rotate=FALSE,useCols=TRUE,points=FALSE){
  
    quartz(width=9,height=9)
    par(mfrow=c(1,1),mar=c(5,5,4,2))
  if(useCols){
    col=coldatamap[toupper(clabels),"color"]
    pnts=pntdatapos[toupper(plabels),"pch"]
  }else{
    col=rainbow(n=length(unique(plabels)))[as.numeric(factor(plabels))]
  }
  
  if(rotate){
    xvals=posdata[as.character(plabels),"Long"]
    yvals=posdata[as.character(plabels),"Lat"]
    phivals=seq(0,2*pi,length=360)

#    calccor=function(p){cor(c(rotPC1(p),rotPC2(p)),c(yvals,xvals))}
#    calccor=function(p){cor(c(s1*rotPC1(p),s2*rotPC2(p)),c(yvals,xvals))}
    calccor=function(p){cor(rotPC1(p,sx,sy),yvals)+cor(rotPC2(p,sx,sy),xvals)}
    corvals=apply(as.matrix(phivals),1,calccor)
    phi=phivals[which.max(corvals)]
 #   phi=-0.279
    
    ang=round(phi/(2*pi)*360)
    ang=min(ang,ang-360)
    
#    plot(rotPC2(phi,sx,sy),rotPC1(phi,sx,sy),type="n",asp=1,frame.plot=FALSE,xlab="\"East-west\" in genetic (PC1-PC2) space",ylab="\"North-south\" in genetic (PC1-PC2) space")
    # par(mar=c(5,2,4,2))
     par(mfrow=c(1,1),mar=c(0,0,0,0),bg="transparent")

     plot(rotPC2(phi,sx,sy),rotPC1(phi,sx,sy),xlab="",ylab="",type="n",asp=1,frame.plot=FALSE,xaxt="n",yaxt="n",xlim=c(-0.07,0.09))

    print(paste("PC",x," rotated by ",ang," degrees.  ","PC",y," rotated by ",ang," degrees (phi=",phi,sep=""))
    x0=0;y0=0;m=0.05;caxis="black" #caxis="darkgrey"

    arrowsOn=FALSE
    if(arrowsOn){
    #NE
    m=0.095
    arrows(0,0,x0+m*cos(-phi),y0+m*sin(-phi),col=caxis,length=0.1)
    #NW
    m=0.07
    arrows(0,0,x0+m*cos(-phi+pi/2),y0+m*sin(-phi+pi/2),col=caxis,length=0.1)
    #SW
    m=0.06 #m=0.065
    arrows(0,0,x0-m*cos(-phi),y0-m*sin(-phi),col=caxis,length=0.1)
    #SE
    m=0.0725
    arrows(0,0,x0-m*cos(-phi+pi/2),y0-m*sin(-phi+pi/2),col=caxis,length=0.1)
  }

#    lines(c(0.075,0.09),c(-0.06,-0.06))
#    lines(c(0.075,0.075),c(-0.0585,-.0615))
#    lines(c(0.09,0.09),c(-0.0585,-.0615))
#    text(c(0.0825),-c(0.0625),c("0.015"),cex=1.25)

    
#    lines(c(-0.05,-0.035),c(0.01,0.01))
#    lines(c(-0.05,-0.05),c(0.0085,.0115))
#    lines(c(-0.035,-0.035),c(0.0085,.0115))
#    text(c(-0.0425),c(0.0125),c("0.015"),cex=1)

    
#    segments(x0-m*cos(-phi),y0-m*sin(-phi),x0+m*cos(-phi),y0+m*sin(-phi),col=caxis)
 #   segments(x0-m*cos(-phi+pi/2),y0-m*sin(-phi+pi/2),x0+m*cos(-phi+pi/2),y0+m*sin(-phi+pi/2),col=caxis)

    
#### OLD VERSION
    pord=sample(length(PCA$ID))
#text(rotPC2(phi,sx,sy)[pord],rotPC1(phi,sx,sy)[pord],alabels[pord],col=col[pord])
    print(phi)
#######

   PC2median=tapply(rotPC2(phi,sx,sy),alabels,median)
   PC1median=tapply(rotPC1(phi,sx,sy),alabels,median)
   colmedian=tapply(col,alabels,function(x){return(x[1])}) # kludge

#   PC2median=tapply(rotPC2(phi,sx,sy)[!OutliersPC1PC2()],alabels[!OutliersPC1PC2()],median)
 #  PC1median=tapply(rotPC1(phi,sx,sy)[!OutliersPC1PC2()],alabels[!OutliersPC1PC2()],median)
 #  colmedian=tapply(col[!OutliersPC1PC2()],alabels[!OutliersPC1PC2()],function(x){return(x[1])}) # kludge

    
  text(rotPC2(phi,sx,sy)[pord],rotPC1(phi,sx,sy)[pord],alabels[pord],col=col[pord],cex=0.75)
#  text(rotPC2(phi,sx,sy)[OutliersPC1PC2(sx,sy)],rotPC1(phi,sx,sy)[OutliersPC1PC2(sx,sy)],alabels[OutliersPC1PC2(sx,sy)],col=1,cex=1)
#  text(rotPC2(phi,sx,sy)[OutliersPC1PC2(sx,sy)],rotPC1(phi,sx,sy)[OutliersPC1PC2(sx,sy)],PCA[OutliersPC1PC2(sx,sy),"ID"],col=1,cex=1)


       # Manual adjustment of specific points
   PC1medianOrig=PC1median
   PC2medianOrig=PC2median
   PC2median[names(PC1median)=="BG"]=PC2median[names(PC1median)=="BG"]-0.0055
   PC2median[names(PC1median)=="MK"]=PC2median[names(PC1median)=="MK"]+0.005
   PC1median[names(PC1median)=="SI"]=PC1median[names(PC1median)=="SI"]-0.005
   PC2median[names(PC1median)=="HR"]=PC2median[names(PC1median)=="HR"]+0.005
   PC1median[names(PC1median)=="RS"]=PC1median[names(PC1median)=="RS"]-0.005
   PC2median[names(PC1median)=="BA"]=PC2median[names(PC1median)=="BA"]+0.005
   PC1median[names(PC1median)=="CY"]=PC1median[names(PC1median)=="CY"]-0.0055
#   PC2median[names(PC2median)=="CY"]=PC2median[names(PC1median)=="CY"]-0.005

   segments(PC2median,PC1median,PC2medianOrig,PC1medianOrig,col=1)
   points(PC2medianOrig,PC1medianOrig,col=colmedian,pch=16,cex=1.5)
   points(PC2median,PC1median,col=colmedian,pch=16,cex=5.5)
   text(PC2median,PC1median,names(PC1median))
  #  recover()

   # Plot legend
   #xleg=rep(0.085,length(PC1median))
   #yleg=seq(-0.07,0.07,length=length(PC1median))
   #points(xleg,yleg,col=colmedian,pch=16,cex=5.5)  
   #text(xleg,yleg,names(PC1median))
   #text(xleg,yleg,levels(plabels))

    
  }else{
    if(points){
     plot(sx*PCA[,paste("PC",x,sep="")],sy*PCA[,paste("PC",y,sep="")],xlab=paste("PC",x,sep=""),ylab=paste("PC",y,sep=""),pch=3,col=col)


     PCxmedian=tapply(PCA[,paste("PC",x,sep="")],plabels,median)
     PCymedian=tapply(PCA[,paste("PC",y,sep="")],plabels,median)
     colmedian=tapply(col,plabels,function(x){return(x[1])}) # kludge
     text(sx*PCxmedian,sy*PCymedian,names(PCxmedian),col=colmedian,)

 
    # recover()
    }else{
    plot(sx*PCA[,paste("PC",x,sep="")],sy*PCA[,paste("PC",y,sep="")],type="n",xlab=paste("PC",x,sep=""),ylab=paste("PC",y,sep=""))
    text(sx*PCA[,paste("PC",x,sep="")],sy*PCA[,paste("PC",y,sep="")],alabels,col=col)
  }
  }

  if(map){
#    par(mar=c(20,3,12,3))
    quartz(height=2.5,width=3)
    par(mar=c(0,0,0,0))
plot(worldpolys,xlim=c(-8,35),ylim=c(35,60),col=coldatamap[toupper(worldpolys$names),"color"])

# tmp=tapply(plabels[alabels!="YG"],factor(alabels[alabels!="YG"]),function(x){return(x[1])})
# colmedian=tapply(col[alabels!="YG"],alabels[alabels!="YG"],function(x){return(x[1])}) # kludge

tmp=tapply(plabels,factor(alabels),function(x){return(x[1])})
colmedian=tapply(col,alabels,function(x){return(x[1])}) # kludge
    

 
 xleg=posdata[tmp,1]
 yleg=posdata[tmp,2]
# recover()
 xleg[tmp=="Russian Federation"]=xleg[tmp=="Russian Federation"]-2
 xleg[tmp=="Albania"]=xleg[tmp=="Albania"]-1
 xleg[tmp=="Slovenia"]=xleg[tmp=="Slovenia"]-1
 yleg[tmp=="Slovenia"]=yleg[tmp=="Slovenia"]-1
 yleg[tmp=="Kosovo"]=yleg[tmp=="Kosovo"]+0.25
 xleg[tmp=="Kosovo"]=xleg[tmp=="Kosovo"]+0.35

 xleg[tmp=="Bosnia and Herzegovina"]=xleg[tmp=="Bosnia and Herzegovina"]-1

 xleg[tmp=="Yugoslavia"]=xleg[tmp=="Serbia and Montenegro"]-1.25
 yleg[tmp=="Yugoslavia"]=yleg[tmp=="Serbia and Montenegro"]
 yleg[tmp=="Serbia and Montenegro"]=yleg[tmp=="Serbia and Montenegro"]+1.25

    points(xleg,yleg,col=colmedian,pch=16,cex=2.95)
 #points(xleg,yleg,col=1,pch=1,cex=5.5)
 text(xleg,yleg,levels(factor(alabels)),cex=0.75)

# text(xleg,yleg,levels(factor(alabels[alabels!="YG"])),cex=0.75)

#colors[countrycolors[match(worldpolys$names,toupper(names(countrycolors)))]],main=paste("PC",x,sep=""))

  }
}

# I don't think this function is relevant to our project. -Peter
plotPCmedianBoot=function(x=1,y=2,sx=1,sy=1,map=TRUE,rotate=FALSE,useCols=TRUE,points=FALSE){
  #pdf("PCboot.pdf",width=6,height=6)
  quartz(width=6,height=6)
  phi=6.001
  plot(rotPC2(phi,sx,sy),rotPC1(phi,sx,sy),type="n",xlab="",ylab="",asp=1,frame.plot=FALSE)

  col=coldatamap[toupper(clabels),"color"]
  colmedian=tapply(col,alabels,function(x){return(x[1])}) # kludge

  x0=0;y0=0;m=0.2;caxis="darkgrey"

  segments(x0-m*cos(-phi),y0-m*sin(-phi),x0+m*cos(-phi),y0+m*sin(-phi),col=caxis)
  segments(x0-m*cos(-phi+pi/2),y0-m*sin(-phi+pi/2),x0+m*cos(-phi+pi/2),y0+m*sin(-phi+pi/2),col=caxis)

   PC2median=tapply(rotPC2(phi,sx,sy),alabels,mean)
   PC1median=tapply(rotPC1(phi,sx,sy),alabels,mean)
   #points(PC2median,PC1median,col=colmedian,pch=16,cex=5.5)
   #text(PC2median,PC1median,names(PC1median))
  
  nboot=10000
  for(i in 1:nboot){

   # Sample
   rotPC1med=tapply(rotPC1(phi,sx,sy),alabels,function(x,i){z=sample(x,replace=TRUE);return(mean(z))}) 
   rotPC2med=tapply(rotPC2(phi,sx,sy),alabels,function(x,i){z=sample(x,replace=TRUE);return(mean(z))}) 


   # Add tests (e.g. counter if Romania > Yugoslav)

   
   # Add points to the graph
   #text(rotPC2med,rotPC1med,i,col=colmedian)
   points(rotPC2med,rotPC1med,cex=0.25,col=colmedian)
 }

  #points(PC2median,PC1median,col=colmedian,pch=16,cex=5.5)
  text(PC2median,PC1median,names(PC1median))

  #dev.off()
}

# See demo.pcvsgeo.R in this repository for a working example of this
# function. -Peter
plotPCdistVsGeoDist=function(sx=1,sy=1,phi=-15){

# PC2median=tapply(rotPC2(phi,sx,sy)[!OutliersPC1PC2()],alabels[!OutliersPC1PC2()],median)
# PC1median=tapply(rotPC1(phi,sx,sy)[!OutliersPC1PC2()],alabels[!OutliersPC1PC2()],median)
# tmp=tapply(plabels[!OutliersPC1PC2()],factor(alabels[!OutliersPC1PC2()]),function(x){return(x[1])}) 
# alabels=alabels[!OutliersPC1PC2()]

 PC2median=tapply(rotPC2(phi,sx,sy),alabels,median)
 PC1median=tapply(rotPC1(phi,sx,sy),alabels,median)
 tmp=tapply(plabels,factor(alabels),function(x){return(x[1])})

 xleg=posdata[tmp,1]
 yleg=posdata[tmp,2]    

 geodistMatrix=apply(cbind(xleg,yleg),1,function(x){return(spDistsN1(cbind(xleg,yleg),x,longlat=TRUE))})
 # Temporarily compute just Euclidean distance (no weighting by eigenvalues)
 PCdistMatrix=apply(cbind(PC1median,PC2median),1,function(x){return(spDistsN1(cbind(PC1median,PC2median),x,longlat=FALSE))})

 # Code for a distance weighted by eigenvalues - slight increase in performance if we scale by inverse of eigenvalue
 weighteddist=function(M,z,l1,l2){
   return(apply(M,1,function(x){return(sqrt(l1*(x[1]-z[1])^2+l2*(x[2]-z[2])^2))}))
 }
 PCdistMatrix=apply(cbind(PC1median,PC2median),1,function(x){return(weighteddist(cbind(PC1median,PC2median),x,l1=1/4.09,l2=1/2.04))})

 rawlabels=expand.grid(levels(factor(alabels)),levels(factor(alabels)))
 pwlabelsMatrix=paste(rawlabels$Var1,rawlabels$Var2)
 n=length(levels(factor(alabels)))
 ids=apply(combn(1:n,2),2,function(x){x[1]+n*(x[2]-1)})

 # Set up new variables 
 geodist=geodistMatrix[ids]
 PCdist=PCdistMatrix[ids]
 pwlabels=pwlabelsMatrix[ids]

 lPC=lm(PCdist~geodist)
 residsPC=abs(lPC$residuals)
 outlier=(residsPC>quantile(residsPC,0.975))
 quartz(width=9,height=9)
 plot(geodist,PCdist,type="n",xlab="Geographic distance (km)",ylab="Euclidean distance in PC1-PC2 space")
 points(geodist[!outlier],PCdist[!outlier],pch=16,cex=1.25)
 text(geodist[outlier],PCdist[outlier],pwlabels[outlier],cex=1)

 alpha=0.025
 abline(lPC)
 abline(lPC$coefficients[1]+quantile(residsPC,1-alpha),lPC$coefficients[2],lty=2)
 abline(lPC$coefficients[1]-quantile(residsPC,1-alpha),lPC$coefficients[2],lty=2)
 print(summary(lPC))

 
 print(cor(geodist,PCdist))
 print(cor(geodist[!outlier],PCdist[!outlier]))
 print(pwlabels[outlier])

}

# This function is marked as "old" to I assume that it is not
# relevant. -Peter
OutliersPC1PC2old=function(sx=-1,sy=-1,phi=-16){

    phi=(phi/360)*2*pi

    xvals=posdata[as.character(plabels),"Long"]
    yvals=posdata[as.character(plabels),"Lat"]

    residsPC1=abs(lm(rotPC1(phi,sx,sy)~yvals)$residuals)
    alpha=0.025
    outliersPC1=(residsPC1>quantile(residsPC1,1-alpha))
    
    residsPC2=abs(lm(rotPC2(phi,sx,sy)~xvals)$residuals)
    alpha=0.025
    outliersPC2=(residsPC2>quantile(residsPC2,1-alpha))

    outliers=outliersPC1|outliersPC2

    outtab=table(factor(alabels[outliers],levels=levels(factor(alabels))))/table(factor(alabels))
    ord=order(outtab,decreasing=T)
    print(cbind(table(factor(alabels[outliers],levels=levels(factor(alabels))))[ord],table(alabels)[ord],signif(outtab[ord],2)))
    
    return(outliers)
}

# It looks like this function is used to "outlying" identify pairs of
# individuals, and does not seem essential to the tools we are
# developing. -Peter
OutliersPC1PC2=function(sx=-1,sy=1,phi=-16,plots=FALSE,alpha=0.025
){

    phi=(phi/360)*2*pi

    xvals0=posdata[as.character(plabels),"Long"]
    yvals0=posdata[as.character(plabels),"Lat"]

   
    xvals=(xvals0-mean(xvals0))/sd(xvals0)
    yvals=(yvals0-mean(yvals0))/sd(yvals0)

    pcNS0=rotPC1(phi,sx,sy)
    pcEW0=rotPC2(phi,sx,sy)

    pcNS=(pcNS0-mean(pcNS0))/sd(pcNS0)
    pcEW=(pcEW0-mean(pcEW0))/sd(pcEW0)

    dvec=sqrt((pcNS-yvals)^2+(pcEW-xvals)^2)

    outliers=dvec>quantile(dvec,1-alpha)


    if(plots){
      #quartz()
      #plot(pcEW,pcNS,type="n")
      #text(pcEW,pcNS,alabels,col=1)
      #text(xvals,yvals,alabels,col=2)

      refpts1=1387  # a UK individual
      refpts2=1195  # an Italina individual
      ScaleFactor=spDistsN1(matrix(c(xvals0[refpts1],yvals0[refpts1]),nrow=1),c(xvals0[refpts2],yvals0[refpts2]),longlat=TRUE)/spDistsN1(matrix(c(xvals[refpts1],yvals[refpts1]),nrow=1),c(xvals[refpts2],yvals[refpts2]))
      print(ScaleFactor)
      quartz()
      hist(dvec,breaks=30,freq=FALSE,xlab="Distance between PC1-PC2 and geographic position",main="")
      abline(v=quantile(dvec,1-alpha),col="red")
      quartz(height=6,width=12)
      split.screen(c(1,2))

      screen(1)
      plot(pcEW0,pcNS0,pch=3,col="grey",xlab="\"East-west\" in rotated PC1-PC2 space",ylab="\"North-south\" in rotated PC1-PC2 space")
      text(pcEW0[outliers],pcNS0[outliers],alabels[outliers],col=2)
      screen(2)
      plot(pcEW0,pcNS0,pch=3,col="grey",xlab="\"East-west\" in rotated PC1-PC2 space",ylab="\"North-south\" in rotated PC1-PC2 space")
      text(pcEW0[outliers],pcNS0[outliers],PCA[outliers,"ID"],col=2)

      quartz()
      dord=order(dvec)
      plot(dvec[dord],seq(0,1,length=length(dvec)),type="n")
      points(dvec[dord][!outliers[dord]],seq(0,1,length=length(dvec))[!outliers[dord]],pch=3,col="gray")
      text(dvec[dord][outliers[dord]],seq(0,1,length=length(dvec))[outliers[dord]],alabels[dord][outliers[dord]])
      abline(v=quantile(dvec,1-alpha),col="red")
    }
 

    
    outtab=table(factor(alabels[outliers],levels=levels(factor(alabels))))/table(factor(alabels))
    ord=order(outtab,decreasing=T)
    print(cbind(table(factor(alabels[outliers],levels=levels(factor(alabels))))[ord],table(alabels)[ord],signif(outtab[ord],2)))
    
    return(outliers)
}

# It looks like this function is used to "outlying" identify pairs of
# individuals, and does not seem essential to the tools we are
# developing. -Peter
OutliersPC1PC2country=function(sx=-1,sy=1,phi=-16,plots=FALSE,alpha=0.025,...
){

    phi=(phi/360)*2*pi

    pcEWmedian0=tapply(rotPC2(phi,sx,sy),alabels,median)
    pcNSmedian0=tapply(rotPC1(phi,sx,sy),alabels,median)
    labels=tapply(plabels,factor(alabels),function(x){return(x[1])})
    labels2=tapply(alabels,factor(alabels),function(x){return(x[1])})


    pcNS=(pcNSmedian0-mean(pcNSmedian0))/sd(pcNSmedian0)
    pcEW=(pcEWmedian0-mean(pcEWmedian0))/sd(pcEWmedian0)


    xvals0=posdata[as.character(labels),"Long"]
    yvals0=posdata[as.character(labels),"Lat"]

    xvals=(xvals0-mean(xvals0))/sd(xvals0)
    yvals=(yvals0-mean(yvals0))/sd(yvals0)

    dvec=sqrt((pcNS-yvals)^2+(pcEW-xvals)^2)

    outliers=dvec>quantile(dvec,1-alpha)

    if(plots){
      refpts1=1  
      refpts2=20  
      ScaleFactor=spDistsN1(matrix(c(xvals0[refpts1],yvals0[refpts1]),nrow=1),c(xvals0[refpts2],yvals0[refpts2]),longlat=TRUE)/spDistsN1(matrix(c(xvals[refpts1],yvals[refpts1]),nrow=1),c(xvals[refpts2],yvals[refpts2]))
      print(ScaleFactor)

      quartz()
      hist(dvec,breaks=10,freq=FALSE,xlab="Distance between PC1-PC2 and geographic position",main="")
      abline(v=quantile(dvec,1-alpha),col="red")


      quartz()
      dord=order(dvec)
      plot(dvec[dord],seq(0,1,length=length(dvec)),type="n",xlab="Distance between PC1-PC2 and geographic position",ylab="Quantile")
      text(dvec[dord],seq(0,1,length=length(dvec)),labels2[dord],cex=0.5)

#      points(dvec[dord][!outliers[dord]],seq(0,1,length=length(dvec))[!outliers[dord]],pch=3,col="gray")
#      text(dvec[dord][outliers[dord]],seq(0,1,length=length(dvec))[outliers[dord]],labels2[dord][outliers[dord]])
      abline(h=1-alpha,col="gray")
      abline(v=quantile(dvec,1-alpha,...),col="red")

      
    }

}

# See demo.pcvslatlong.R in this repository for a working example of
# this function. -Peter
plotPCvsLatLong=function(sx=-1,sy=-1,phi=-16){

    phi=(phi/360)*2*pi
    print(phi)
    xvals=posdata[as.character(plabels),"Long"]
    yvals=posdata[as.character(plabels),"Lat"]
    col=coldatamap[toupper(clabels),"color"]
    
    
    quartz(height=5,width=10)
    par(mfrow=c(1,2))

    plot(yvals,rotPC1(phi,sx,sy),xlab="Latitude",ylab="\"North-south\" in rotated PC1-PC2 space",pch=16,type="n")
    lPC1=lm(rotPC1(phi,sx,sy)~yvals)
    residsPC1=abs(lPC1$residuals)
    alpha=0.025
    outliers=OutliersPC1PC2(sx,sy);#(residsPC1>quantile(residsPC1,1-alpha))
#    points(yvals[!outliers],sy*rotPC1(phi)[!outliers],pch=16,col=col[!outliers])
    points(yvals,rotPC1(phi,sx,sy),pch=16,col=col)
    text(yvals[outliers],rotPC1(phi,sx,sy)[outliers],alabels[outliers],cex=0.75)

    abline(lPC1)
    abline(lPC1$coefficients[1]+quantile(residsPC1,1-alpha),lPC1$coefficients[2],lty=2)
    abline(lPC1$coefficients[1]-quantile(residsPC1,1-alpha),lPC1$coefficients[2],lty=2)
    print(summary(lPC1))
    print(paste("Lat Correlation = ",signif(cor(yvals,rotPC1(phi,sx,sy)),3)))
    print(paste("Lat Correlation (no outliers) = ",signif(cor(yvals[!outliers],rotPC1(phi,sx,sy)[!outliers]),3)))

    
    plot(xvals,rotPC2(phi,sx,sy),xlab="Longitude",ylab="\"East-west\" in rotated PC1-PC2 space",type="n")
    lPC2=lm(rotPC2(phi,sx,sy)~xvals)
    residsPC2=abs(lPC2$residuals)
    alpha=0.025
    outliers=OutliersPC1PC2(sx,sy);#(residsPC2>quantile(residsPC2,1-alpha))
    points(xvals,rotPC2(phi,sx,sy),pch=16,col=col)
    text(xvals[outliers],rotPC2(phi,sx,sy)[outliers],alabels[outliers],cex=0.75)

    abline(lPC2)
    abline(lPC2$coefficients[1]+quantile(residsPC2,1-alpha),lPC2$coefficients[2],lty=2)
    abline(lPC2$coefficients[1]-quantile(residsPC2,1-alpha),lPC2$coefficients[2],lty=2)
    print(summary(lPC1))

    print(summary(lPC2))
    print(paste("Long Correlation = ",signif(cor(xvals,rotPC2(phi,sx,sy)),3)))
    print(paste("Long Correlation (no outliers) = ",signif(cor(xvals[!outliers],rotPC2(phi,sx,sy)[!outliers]),3)))

    print(summary(lm(PCA[,"PC2"]~xvals)))
    print(summary(lm(PCA[,"PC1"]~yvals)))

}

# See demo.pcvslatlongmedian.R in this repository for a working
# example of this function. -Peter
plotPCvsLatLongMedian=function(sx=-1,sy=-1,phi=-16){

    phi=(phi/360)*2*pi
    print(phi)

    PC2median=tapply(rotPC2(phi,sx,sy),maplabel,median)
    PC1median=tapply(rotPC1(phi,sx,sy),maplabel,median)
    PC2se=tapply(rotPC2(phi,sx,sy),maplabel,function(x){sd(x)/sqrt(length(x))})
    PC1se=tapply(rotPC1(phi,sx,sy),maplabel,function(x){sd(x)/sqrt(length(x))})
    labels=tapply(plabels,factor(maplabel),function(x){return(x[1])})

    xvals=posdata[as.character(levels(factor(maplabel))),"Long"]
    yvals=posdata[as.character(levels(factor(maplabel))),"Lat"]

  #  xvals=coordinates(worldpolys)[match(toupper(levels(factor(maplabel))),worldpolys$names),1]
  #  yvals=coordinates(worldpolys)[match(toupper(levels(factor(maplabel))),worldpolys$names),2]

    posdata=data.frame(xvals,yvals)
    names(posdata)=c("Long","Lat")
    row.names(posdata)=levels(factor(maplabel))
    print(posdata)

    tmp=tapply(plabels,factor(maplabel),function(x){return(x[1])}) 
    xleg=posdata[tmp,1]
    yleg=posdata[tmp,2]    


    col=coldatamap[toupper(clabels),"color"]
    colmedian=tapply(col,maplabel,function(x){return(x[1])}) # kludge
    labelsA=tapply(alabels,maplabel,function(x){return(x[1])}) # kludge

#    PC2median=tapply(rotPC2(phi,sx,sy)[!OutliersPC1PC2()],alabels[!OutliersPC1PC2()],median)
#    PC1median=tapply(rotPC1(phi,sx,sy)[!OutliersPC1PC2()],alabels[!OutliersPC1PC2()],median)
#    tmp=tapply(plabels[!OutliersPC1PC2()],factor(alabels[!OutliersPC1PC2()]),function(x){return(x[1])}) 

    print(paste(length(xleg),length(yleg),length(bboxes$S),length(PC1median)))
    print(PC1median)
    quartz(height=5,width=10)
    par(mfrow=c(1,2))

    plot(yleg,PC1median,xlab="Latitude",ylab="\"South-North\" in rotated PC1-PC2 space",pch=16,type="n")

    lPC1=lm(PC1median~yleg)
    abline(lPC1)
    
    residsPC1=abs(lPC1$residuals)
    apply(cbind(yleg,PC1median,2*PC1se),1,function(x){lines(c(x[1],x[1]),c(x[2]-x[3],x[2]+x[3]));lines(c(x[1]-0.5,x[1]+0.5),c(x[2]+x[3],x[2]+x[3]));lines(c(x[1]-0.5,x[1]+0.5),c(x[2]-x[3],x[2]-x[3]))})
    apply(cbind(PC1median,bboxes$S,bboxes$N),1,function(x){lines(c(x[2],x[3]),c(x[1],x[1]));lines(c(x[2],x[2]),c(x[1]-0.0025,x[1]+0.0025));lines(c(x[3],x[3]),c(x[1]-0.0025,x[1]+0.0025))})

    points(yleg,PC1median,pch=16,cex=1.5,col=colmedian)
    text(yleg,PC1median,labelsA,cex=0.75)

    alpha=0.025
  #  outliers=residsPC1>quantile(residsPC1,1-alpha)
  #  points(yleg[!outliers],PC1median[!outliers],pch=16,cex=1.5,col=colmedian[!outliers]) 
  #  text(yleg[outliers],PC1median[outliers],labelsA[outliers])
    abline(lPC1)
    abline(lPC1$coefficients[1]+quantile(residsPC1,1-alpha),lPC1$coefficients[2],lty=2)
    abline(lPC1$coefficients[1]-quantile(residsPC1,1-alpha),lPC1$coefficients[2],lty=2)

#    abline(lPC1$coefficients[1]+2*sqrt(median(residsPC1)^2),lPC1$coefficients[2],lty=2,col=2)
#    abline(lPC1$coefficients[1]-2*sqrt(median(residsPC1)^2),lPC1$coefficients[2],lty=2,col=2)

    print(summary(lPC1))
    
 #   text(yvals[outliers],sy*rotPC1(phi)[outliers],alabels[outliers])


    plot(xleg,PC2median,xlab="Longitude",ylab="\"West-East\" in rotated PC1-PC2 space",pch=16,type="n")
    lPC2=lm(PC2median~xleg)
    residsPC2=abs(lPC2$residuals)
    
    apply(cbind(xleg,PC2median,2*PC2se),1,function(x){lines(c(x[1],x[1]),c(x[2]-x[3],x[2]+x[3]));lines(c(x[1]-0.5,x[1]+0.5),c(x[2]+x[3],x[2]+x[3]));lines(c(x[1]-0.5,x[1]+0.5),c(x[2]-x[3],x[2]-x[3]))})
    apply(cbind(PC2median,bboxes$W,bboxes$E),1,function(x){lines(c(x[2],x[3]),c(x[1],x[1]));lines(c(x[2],x[2]),c(x[1]-0.0025,x[1]+0.0025));lines(c(x[3],x[3]),c(x[1]-0.0025,x[1]+0.0025))})

    points(xleg,PC2median,pch=16,cex=1.5,col=colmedian)
    text(xleg,PC2median,labelsA,cex=0.75)

#    points(xleg[!outliers],PC2median[!outliers],pch=16)
#    text(xleg[outliers],PC2median[outliers],names(PC2median)[outliers])
#    text(xleg[names(PC1median)=="RO"],PC2median[names(PC1median)=="RO"],"RO")
#    text(xleg[names(PC1median)=="YG"],PC2median[names(PC1median)=="YG"],"YG")

    abline(lPC2)
    abline(lPC2$coefficients[1]+quantile(residsPC2,1-alpha),lPC2$coefficients[2],lty=2)
    abline(lPC2$coefficients[1]-quantile(residsPC2,1-alpha),lPC2$coefficients[2],lty=2)

#    abline(lPC2$coefficients[1]+2*sqrt(median(residsPC2)^2),lPC2$coefficients[2],lty=2,col=2)
#    abline(lPC2$coefficients[1]-2*sqrt(median(residsPC2)^2),lPC2$coefficients[2],lty=2,col=2)

    print(summary(lPC2))

}


# This seems like a specialized function and not within the scope of
# the package we are developing. -Peter
plotPCvarVsGeoSpread=function(sx=-1,sy=-1,phi=-16){

  NSvar=tapply(rotPC1(phi,sx,sy),maplabel,var)
  EWvar=tapply(rotPC2(phi,sx,sy),maplabel,var)

  xvals=posdata[as.character(levels(factor(maplabel))),"Long"]
  yvals=posdata[as.character(levels(factor(maplabel))),"Lat"]
  tmp=tapply(plabels,factor(maplabel),function(x){return(x[1])}) 
  xleg=posdata[tmp,1]
  yleg=posdata[tmp,2]    

  col=coldatamap[toupper(clabels),"color"]
  colmean=tapply(col,maplabel,function(x){return(x[1])}) # kludge
  labelsA=tapply(alabels,maplabel,function(x){return(x[1])}) # kludge

  NSdist=apply(matrix(levels(factor(maplabel)),ncol=1),1,function(x){spDistsN1(matrix(c(posdata[as.character(x),"Long"],bboxes[as.character(x),"S"]),nrow=1),c(posdata[as.character(x),"Long"],bboxes[as.character(x),"N"]),longlat=TRUE)})
  EWdist=apply(matrix(levels(factor(maplabel)),ncol=1),1,function(x){spDistsN1(matrix(c(posdata[as.character(x),"Lat"],bboxes[as.character(x),"E"]),nrow=1),c(posdata[as.character(x),"Lat"],bboxes[as.character(x),"W"]),longlat=TRUE)})
  print(cbind(levels(factor(maplabel)),NSdist))
  print(cbind(levels(factor(maplabel)),EWdist))

  quartz()
  par(mfrow=c(3,1))

  samsizeok=(table(maplabel)>10)#&(labelsA!="IT")
  print(sum(samsizeok))
                                        #  plot(NSdist[samsizeok],log(NSvar[samsizeok]),xlim=c(0,1750),ylim=c(-11.5,log(3e-4)),type="n",xlab="",ylab="")
#  text(NSdist[samsizeok],log(NSvar[samsizeok]),labelsA[samsizeok],col=2)
#  text(EWdist[samsizeok],log(EWvar[samsizeok]),labelsA[samsizeok],col=3)
#  abline(lm(log(NSvar[samsizeok])~NSdist[samsizeok]))
#  abline(lm(log(EWvar[samsizeok])~EWdist[samsizeok]))
#  title(xlab="Distance (km)",ylab="Variation in PC coordinates")

  hist(log10(EWvar[samsizeok]))
  hist(log10(EWdist[samsizeok]))
  
  plot(log10(NSdist[samsizeok]),log10(NSvar[samsizeok]),xlim=c(log10(200),log10(1750)),ylim=c(log10(0.1e-4),log10(3e-4)),type="n",xlab="",ylab="")
  text(log10(NSdist[samsizeok]),log10(NSvar[samsizeok]),labelsA[samsizeok],col="red")
  text(log10(EWdist[samsizeok]),log10(EWvar[samsizeok]),labelsA[samsizeok],col="blue")
  abline(lm(log10(NSvar[samsizeok])~log10(NSdist[samsizeok])),col="red")
  abline(lm(log10(EWvar[samsizeok])~log10(EWdist[samsizeok])),col="blue")
  title(xlab="log10[Distance (km)]",ylab="Variation in PC coordinates")

  print(summary(lm(NSvar[samsizeok]~NSdist[samsizeok])))
  print(summary(lm(EWvar[samsizeok]~EWdist[samsizeok])))
  print("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX")
#  print(summary(lm(log10(NSvar[samsizeok])~(NSdist[samsizeok]))))
#  print(summary(lm(log10(EWvar[samsizeok])~(EWdist[samsizeok]))))
  print(summary(lm(log10(NSvar[samsizeok])~log10(NSdist[samsizeok]))))
  print(summary(lm(log10(EWvar[samsizeok])~log10(EWdist[samsizeok]))))
  print(cor.test(log10(EWvar[samsizeok]),log10(EWdist[samsizeok]),method="spearman"))
  print(cor.test(log10(NSvar[samsizeok]),log10(NSdist[samsizeok]),method="spearman"))
#  plot(log(NSvar[samsizeok]/EWvar[samsizeok])~log(NSdist[samsizeok]/EWdist[samsizeok]))
  print(summary(lm(log(NSvar[samsizeok]/EWvar[samsizeok])~log(NSdist[samsizeok]/EWdist[samsizeok]))))
  quartz()
  plot(lm(log10(EWvar[samsizeok])~log10(EWdist[samsizeok])))
  
#  print(summary(lm(log(NSvar)~(NSdist))))
#  print(summary(lm(log(EWvar)~(EWdist))))
#  print(summary(lm(log(NSvar[samsizeok])~(NSdist[samsizeok]))))
#  print(summary(lm(log(EWvar[samsizeok])~(EWdist[samsizeok]))))

                                        #  print(summary(lm((NSvar[samsizeok])~log(NSdist[samsizeok]))))
#  print(summary(lm((EWvar[samsizeok])~log(EWdist[samsizeok]))))

  
#  EWdist=apply(maplabel,1,fun
  
#    apply(cbind(PC1mean,bboxes$S,bboxes$N),1,function(x){lines(c(x[2],x[3]),c(x[1],x[1]));lines(c(x[2],x[2]),c(x[1]-0.0025,x[1]+0.0025));lines(c(x[3],x[3]),c(x[1]-0.0025,x[1]+0.0025))})

}

# See demo.plotpc.R in this repository for a working example of this
# function. -Peter
plotPC=function(x,sx=1){

aboveThresh=levels(factor(clabels))[table(clabels)>=2]
labelsThin=clabels[is.element(clabels,aboveThresh)]

PCvals=PCA[is.element(clabels,aboveThresh),paste("PC",x,sep="")]*sx
meanvals=tapply(PCvals,factor(labelsThin),median)
sdvals=tapply(PCvals,factor(labelsThin),sd)

print(summary(PCvals))
bins=seq(1.01*min(PCvals),1.01*max(PCvals),length=15)
#bins=seq(0.99*min(PCvals),1.01*max(PCvals),length=18)
#bins=seq(-0.085,0.085,length=20)
colors=rainbow(length(bins)-1,s=0.6,start=0/6,end=4/6)
colors=heat.colors(length(bins)-1)


countrycolors=apply(meanvals,1,function(x){for(i in 1:length(bins)) if(x<bins[i]){return(i-1); break;}})


par(fin=c(10,8))
layout(matrix(c(1,2,3),1,3),widths=c(0.75,0.05,0.2))

par(mar=c(5,4,4,2))
plot(worldpolys,xlim=c(-12,55),ylim=c(35,65),col=colors[countrycolors[match(worldpolys$names,toupper(names(countrycolors)))]],main=paste("PC",x,sep=""))
xvals=coordinates(worldpolys)[match(toupper(names(countrycolors)),worldpolys$names),1]
yvals=coordinates(worldpolys)[match(toupper(names(countrycolors)),worldpolys$names),2]
#recover()
#text(xvals,yvals,paste(signif(2*sdvals,digits=1)),col=2,cex=cex.text*1.5)

#paste(signif(meanvals,digits=1),"\n(",signif(2*sdvals,digits=1), ")")

table(factor(clabels))
par(mar=c(8,0.5,8,0.25))
image(y=bins[1:(length(bins)-1)], z=matrix(bins[1:(length(bins)-1)], ncol=length(bins)-1), col=colors, axes=FALSE, ylab="", xlab="")
#axis(4,at=signif(bins,digits=2))
binmarks=signif(bins,digits=2)
binlabels=paste("[ ",binmarks[1:length(binmarks)]," , ",binmarks[2:length(binmarks)],")",sep="")
axis(4,at=bins,labels=binlabels,las=2,cex.axis=cex.text*1.5)
#level <- (bins[2]-bins[1])/2
#segments(-1, level, 1, level, col="black", lwd=1)
box()

# plot the axis for the significance bar
#par(mar=c(10.95, 0, 10.95, 0))
#plot(0:1, range(bins), axes=F, ylim=range(bins),xlab="", ylab="", type="n")
#segments(0, bins, 0.08, bins)
#text(0.2, bins,
#   labels=paste(c(rep("", max(bins))), signif(bins,digits=2),
#   sep=""), adj=c(0,0.5), cex=cex.text*1.5)

#par(mar=rep(0,4))
#plot(0:1, 0:1, xlab="", ylab="", type="n", axes=F)
#text(0.25, 1, labels="p-value", adj=0.5, cex=cex.text)

}

# Unclear what this function does; difficult to reproduce steps. -Peter
plotPCstruct=function(x,colopt=1,sx=1,ordopt=1,spacing=1){

labels=rlabels

#zchar=demodata[as.character(PCA$ID),"Wave"]
#zchar=demodata[as.character(PCA$ID),"tagGPnodiscord"]==demodata[as.character(PCA$ID),"Country_Self"]
#ddEuro=demodata[as.character(PCA$ID),]
#zchar=!(ddEuro$Country_Self==ddEuro$tagGPnodiscord | (is.element(ddEuro$tagGPnodiscord,c("Swiss-French","Swiss-German","Swiss-Italian"))&ddEuro$Country_Self=="Switzerland"))

if(ordopt==1){
zchar=PCA[,paste("PC",x,sep="")]*sx
#zchar=ddata[as.character(PCA$ID),"groups3"]

ord=order(grplabels,as.character(labels),zchar)

pops=levels(factor(labels))
grpIDs=factor(as.character(grps[as.character(pops),"grp"]),levels=c("EuropeNW","EuropeN","EuropeNE","EuropeW","EuropeC","EuropeE","EuropeSW","EuropeS","EuropeSE","EuropeESE"),ordered=T)
popsord=pops[order(grpIDs,pops)]

}else if(ordopt==2){

zchar=PCA[,paste("PC",x,sep="")]*sx
pops=levels(factor(labels))
meanvals=tapply(zchar,factor(labels),mean)
popsord=pops[order(meanvals)]
ord=order(meanvals[as.character(labels)],zchar)

} else if(ordopt==3){
  zchar=PCA[,paste("PC",x,sep="")]*sx
  pops=levels(factor(labels))
  xvals=posdata[as.character(levels(factor(plabels))),"Long"]
  popsord=pops[order(xvals)]
  ord=order(posdata[as.character(plabels),"Long"],zchar)
 
}

if(colopt==1){
cols=brewer.pal("Paired",n=7)[as.numeric(factor(demodata[as.character(PCA$ID),"Wave"][ord]))]

}else{
#cols=brewer.pal("Paired",n=10)[match(labels[ord],popsord)%%10+1]
cols=brewer.pal("Paired",n=3)[match(labels[ord],popsord)%%2+1]

}

par(mar=c(10,5,6,2))
n=length(labels)
breaks=c(0,which(labels[ord][1:(n-1)]!=labels[ord][2:n]),n)
nbrks=length(breaks)
midpts=(breaks[1:(nbrks-1)]+breaks[2:nbrks])/2

if(ordopt==1){
gbreaks=c(0,which(grplabels[ord][1:(n-1)]!=grplabels[ord][2:n]),n)
gnbrks=length(gbreaks)
gmidpts=(gbreaks[1:(gnbrks-1)]+gbreaks[2:gnbrks])/2
}

par(mfrow=c(1,1))
plot(1:n,PCA[ord,paste("PC",x,sep="")]*sx,xlim=c(n/25,n-n/25),type="n",xaxt="n",xlab="",ylab=paste("PC",x,sep=""),xpd=FALSE)

#cols=brewer.pal("Paired",n=7)[as.numeric(factor(zchar[ord]))]
#cols=coldata[toupper(clabels[ord]),"color"]
#cols=brewer.pal("Paired",n=7)[as.numeric(cut(zchar[ord],8))]
#cols=heat.colors(n=10)[as.numeric(cut(zchar[ord],10))]
#cols=brewer.pal("Paired",n=10)[match(labels[ord],popsord)%%10+1]
if(spacing==1){
segments(breaks,rep(-1,n),breaks,rep(1,n),lwd=0.5,col="gray")
  rect(0:(n-1),rep(0,n),1:n,PCA[ord,paste("PC",x,sep="")]*sx,col=cols,border=cols)
mtext(popsord,side=1,at=midpts,las=2)
}else{
  cnts=(breaks[2:nbrks]-breaks[1:(nbrks-1)])
  sc=1/cnts  # scale for each population
  scbar=rep(sc,times=cnts)  # scale for each bar
  scbar=(scbar/sum(scbar))*n # renormalize to n
  segments(seq(0,n,length=nbrks),rep(-1,n),seq(0,n,length=nbrks),rep(1,n),lwd=0.5,col="lightgray")
  rect(c(0,cumsum(scbar[1:(n-1)])),rep(0,n),cumsum(scbar),PCA[ord,paste("PC",x,sep="")]*sx,col=cols,border=cols)
  mtext(popsord,side=1,at=seq(0,n,length=nbrks)[1:(nbrks-1)]+(n/nbrks)/2,las=2)

}

if(ordopt==1 & spacing==1){
axis(side=3,at=gbreaks,labels=FALSE,lwd=2,col=2)
mtext(levels(grpIDs),side=3,at=gmidpts,las=2)
}
#recover()
}

# See demo.swiss.R in this repository for a working example of this
# function. -Peter
plotSwiss = function(sx=1,sy=-1,phi=-16) {

quartz(width=6,height=6)
plot(rotPC2(phi,sx,sy)[!is.element(alabels,c("CH-F","CH-I","CH-G"))],rotPC1(phi,sx,sy)[!is.element(alabels,c("CH-F","CH-I","CH-G"))],xlab="",ylab="",asp=1,pch=16,cex=0.5,col="lightgrey",xlim=c(-0.03,0.03),ylim=c(-0.03,0.03))

points(rotPC2(phi,sx,sy)[alabels=="IT"],rotPC1(phi,sx,sy)[alabels=="IT"],col="lightblue",cex=1,pch=17)
points(rotPC2(phi,sx,sy)[alabels=="DE"],rotPC1(phi,sx,sy)[alabels=="DE"],col="pink",cex=1.2,pch=18)
points(rotPC2(phi,sx,sy)[alabels=="FR"],rotPC1(phi,sx,sy)[alabels=="FR"],col="greenyellow",cex=1,pch=15)

points(rotPC2(phi,sx,sy)[alabels=="CH-I"],rotPC1(phi,sx,sy)[alabels=="CH-I"],col="blue",cex=1.5,pch=2)
points(rotPC2(phi,sx,sy)[alabels=="CH-G"],rotPC1(phi,sx,sy)[alabels=="CH-G"],col="red",cex=1.5,pch=5)
points(rotPC2(phi,sx,sy)[alabels=="CH-F"],rotPC1(phi,sx,sy)[alabels=="CH-F"],col="green4",cex=1.5,pch=22)

}

# I don't think this function is relevant to our project. -Peter
anovaAnalysis=function(){
# Analysis of variance on Genotyping wave
plabel=as.factor(demodata[as.character(PCA$ID), "tagGPnodiscord"])
pPC1means=tapply(PCA$PC1,plabel,mean)
pPC2means=tapply(PCA$PC2,plabel,mean)
wave=as.factor(demodata[as.character(PCA$ID), "Wave"])
src=as.factor(demodata[as.character(PCA$ID), "Source"])
m1=PCA$PC1-pPC1means[plabel]~wave
m2=PCA$PC2-pPC2means[plabel]~wave
plot(m1)
plot(m2)
summary(lm(m1))
summary(lm(m2))

# Alternative approach
anova(lm(PCA$PC1~plabel+wave))
anova(lm(PCA$PC2~plabel+wave))
anova(lm(PCA$PC1~plabel+src))
anova(lm(PCA$PC2~plabel+src))
anova(lm(PCA$PC1~plabel+src+wave))
anova(lm(PCA$PC2~plabel+src+wave))

}
