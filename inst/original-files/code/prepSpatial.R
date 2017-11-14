


library(MASS)
library(maps)
library(sp)
library(rgdal)




if(!exists("world_countries")){
  load("world_countries.rda")
}

if(!exists("ddata")){
  ddata=read.table("PCA.txt",header=T,sep="\t",as.is=TRUE)
  row.names(ddata)=ddata$POPRES_ID
  #demodata=read.table("../demodata/d3demodata_V4.txt",header=T,sep="\t",as.is=TRUE)
  #row.names(demodata)=demodata$POPRES_ID

}



fixlabelsA=function(x){
  # Fixed labels for using on maps
x[grep("Wales",x)]="United Kingdom"
x[grep("Scotland",x)]="United Kingdom"
x[grep("Yugoslavia",x)]="Serbia and Montenegro"
x[grep("Serbia",x)]="Serbia and Montenegro"
x[grep("Bosnia",x)]="Bosnia and Herzegovina"
x[grep("Bosnia-Herzegovina",x)]="Bosnia and Herzegovina"
x[grep("Kosovo",x)]="Serbia and Montenegro"
x[grep("Russia",x)]="Russian Federation"
x[grep("Macedonia",x)]="Macedonia, The Former Yugoslav Republic Of"
x[grep("Swiss-German",x)]="Switzerland"
x[grep("Swiss-Italian",x)]="Switzerland"
x[grep("Swiss-French",x)]="Switzerland"
return(x)
}

fixlabelsB=function(x){
  # Fixes labels to get unique geographic positions
x[grep("Yugoslavia",x)]="Serbia and Montenegro"
x[grep("Serbia",x)]="Serbia and Montenegro"
x[grep("Bosnia",x)]="Bosnia and Herzegovina"
x[grep("Bosnia-Herzegovina",x)]="Bosnia and Herzegovina"
#x[grep("Kosovo",x)]="Serbia and Montenegro"
x[grep("Russia",x)]="Russian Federation"
x[grep("Macedonia",x)]="Macedonia, The Former Yugoslav Republic Of"
return(x)
}

abbrvlabels=function(x){
x[grep("United Kingdom",x)]="GB"
x[grep("Wales",x)]="Wales"
x[grep("Scotland",x)]="Sct"
x[grep("Czech Republic",x)]="CZ"
x[grep("Yugoslavia",x)]="YG"
x[grep("Bosnia-Herzegovina",x)]="BA"
x[grep("Kosovo",x)]="KS"
x[grep("Russia",x)]="RU"
x[grep("Iran",x)]="IR"
x[grep("Macedonia",x)]="MK"
x[grep("Swiss-German",x)]="CH"
x[grep("Swiss-Italian",x)]="CH"
x[grep("Swiss-French",x)]="CH"
x[grep("Switzerland",x)]="CH"
x[grep("Spain",x)]="ES"
x[grep("Italy",x)]="IT"
x[grep("Portugal",x)]="PT"
x[grep("France",x)]="FR"
x[grep("Germany",x)]="DE"
x[grep("Greece",x)]="GR"
x[grep("Ireland",x)]="IE"
x[grep("Romania",x)]="RO"
x[grep("Austria",x)]="AT"
x[grep("Hungary",x)]="HU"
x[grep("Belgium",x)]="BE"
x[grep("Poland",x)]="PL"
x[grep("Netherlands",x)]="NL"
x[grep("Norway",x)]="NO"
x[grep("Sweden",x)]="SE"
x[grep("Finland",x)]="FI"
x[grep("Latvia",x)]="LV"
x[grep("Turkey",x)]="TR"
x[grep("Croatia",x)]="HR"
x[grep("Albania",x)]="AL"
x[grep("Ukraine",x)]="UA"
x[grep("Bulgaria",x)]="BG"
x[grep("Slovenia",x)]="SI"
x[grep("Slovakia",x)]="SK"
x[grep("Cyprus",x)]="CY"
x[grep("Denmark",x)]="DK"
x[grep("Serbia",x)]="RS"


return(x)
}

# Different labels

# Labels for mapping : group into polygons
maplabel=fixlabelsA(ddata$tagGPnodiscord)
  
# labels for position-based analyes: based on geographic points
poslabel=fixlabelsB(ddata$tagGPnodiscord)

# Labels for coloring assignment plots
collabelpos=poslabel

collabelmap=poslabel

abbrvlabel=abbrvlabels(ddata$tagGPnodiscord)

############################

# Set-up position data
xvals=coordinates(world_countries)[match(toupper(levels(factor(poslabel))),world_countries$names),1]
yvals=coordinates(world_countries)[match(toupper(levels(factor(poslabel))),world_countries$names),2]


posdata=data.frame(xvals,yvals)
names(posdata)=c("Long","Lat")
row.names(posdata)=levels(factor(poslabel))

if(!exists("bboxes")){
bboxes=apply(matrix(1:length(levels(factor(maplabel)))),1,function(x){bbox(world_countries[match(toupper(levels(factor(maplabel))[x]),world_countries$names),])})
bboxes=data.frame(t(bboxes))
row.names(bboxes)=levels(factor(maplabel))
names(bboxes)=c("W","S","E","N")

bboxes["Norway","W"]=4.75
bboxes["Norway","N"]=71.25
bboxes["United Kingdom","N"]=59
bboxes["United Kingdom","W"]=-8.5
bboxes["United Kingdom","S"]=50
bboxes["Spain","W"]=-9.25
bboxes["Spain","S"]=36
bboxes["Spain","E"]=3.5
bboxes["Portugal","S"]=37
bboxes["Portugal","W"]=-9.25
bboxes["France","W"]=-4.75
bboxes["France","S"]=42.25
bboxes["France","E"]=8
bboxes["Italy","S"]=36.
bboxes["Russian Federation","W"]=28
bboxes["Russian Federation","E"]=55
bboxes["Russian Federation","N"]=70
  

#bboxes["Norway",]=c("NA","NA","NA","NA")
#bboxes["Sweden",]=c("NA","NA","NA","NA")
#bboxes["Finland",]=c("NA","NA","NA","NA")
#bboxes["Slovakia",]=c("NA","NA","NA","NA")
#bboxes["Ukraine",]=c("NA","NA","NA","NA")
#bboxes["Latvia",]=c("NA","NA","NA","NA")
#bboxes["Slovenia",]=c("NA","NA","NA","NA")
#bboxes["Kosovo",]=c("NA","NA","NA","NA")
#bboxes["Bulgaria",]=c("NA","NA","NA","NA")
#bboxes["Albania",]=c("NA","NA","NA","NA")
quartz()
plot(world_countries,xlim=c(-12,55),ylim=c(30,75),col=coldatamap[toupper(world_countries$names),"color"])
#    abline(h=seq(35,65,by=3),col=2)
#    abline(v=seq(-12,55,by=3),col=2)
    apply(bboxes,1,function(x){lines(c(x[1],x[1],x[3],x[3],x[1]),c(x[2],x[4],x[4],x[2],x[2]))})

}


# Corrections
posdata["Scotland","Lat"]=55.95;posdata["Scotland","Long"]=-3.2;
posdata["Switzerland","Lat"]=46.2;posdata["Switzerland","Long"]=6.15; # Note: temp fix for using an old file 
posdata["Swiss-French","Lat"]=46.2;posdata["Swiss-French","Long"]=6.15;
posdata["Swiss-German","Lat"]=47.36;posdata["Swiss-German","Long"]=8.55;
posdata["Swiss-Italian","Lat"]=46;posdata["Swiss-Italian","Long"]=8.95;
posdata["Kosovo","Lat"]=42.7;posdata["Kosovo","Long"]=21.1; # Pristina
posdata["Russian Federation","Lat"]=55.75;posdata["Russian Federation","Long"]=37.5
posdata["Sweden","Lat"]=59.35;posdata["Sweden","Long"]=18 
posdata["Norway","Lat"]=59.93;posdata["Norway","Long"]=10.68 
posdata["Finland","Lat"]=60.17;posdata["Finland","Long"]=24.93 # Helsinki

posdata["Croatia","Lat"]=45.32;posdata["Croatia","Long"]=16.1 # Central
posdata["Italy","Lat"]=42;posdata["Italy","Long"]=12.5  # Rome

#########################

# Set-up color data
coldatamap=read.table("../demodata/ColorTablePCmap.txt",sep="\t",as.is=TRUE)
names(coldatamap)=c("cntry","color")
row.names(coldatamap)=toupper(coldatamap$cntry)

coldatapos=read.table("../demodata/ColorTableAssign.txt",sep="\t",as.is=TRUE)
names(coldatapos)=c("cntry","color")
row.names(coldatapos)=coldatapos$cntry

pntdatapos=read.table("../demodata/PointsTablePCmap.txt",sep="\t",as.is=TRUE)
names(pntdatapos)=c("cntry","pch")
row.names(pntdatapos)=toupper(pntdatapos$cntry)



########################

# Set up order of individuals
grps=read.table("../demodata/groupingsD3_JN3.txt",sep="\t",header=T);names(grps)=c("ID","grp");row.names(grps)=grps$ID

# set up order of groups
grplabels=factor(as.character(ddata[,"groups3"]),levels=c("EuropeNW","EuropeN","EuropeNE","EuropeW","EuropeC","EuropeE","EuropeSW","EuropeS","EuropeSE","EuropeESE"),ordered=T)

########################


plotGeoPoints=function(){
plot(world_countries,xlim=c(-12,50),ylim=c(35,65))
points(posdata$Long,posdata$Lat,pch=16,col=coldatapos[row.names(posdata),"color"])
}

plotGeoMaps=function(){

}

rotPC1=function(p,sx=1,sy=1,x=0,y=0){
   cos(p)*(sy*PCA[,"PC1"]-y)-sin(p)*(sx*PCA[,"PC2"]-x)      
}
rotPC2=function(p,sx=1,sy=1,x=0,y=0){
  cos(p)*(sx*PCA[,"PC2"]-x)+sin(p)*(sy*PCA[,"PC1"]-y)
}
