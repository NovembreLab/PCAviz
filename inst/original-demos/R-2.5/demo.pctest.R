library(sp)
library(rgdal)
load("../../original-files/data/world_countries.rda")
worldpolys=as.data.frame(world_countries)
coldatamap=read.table("../../original-files/data/ColorTablePCmap.txt",sep="\t",as.is=TRUE)
names(coldatamap)=c("cntry","color")
row.names(coldatamap)=toupper(coldatamap$cntry)
plot(worldpolys,xlim=c(-8,35),ylim=c(35,60),col=coldatamap[toupper(worldpolys$names),"color"])

