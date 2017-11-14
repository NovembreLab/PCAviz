# This script only plots a world map; not particularly interesting. -Peter
load("world_countries.rda")
coldatamap=read.table("ColorTablePCmap.txt",sep="\t",as.is=TRUE)
names(coldatamap)=c("cntry","color")
row.names(coldatamap)=toupper(coldatamap$cntry)
plot(worldpolys,xlim=c(-8,35),ylim=c(35,60),col=coldatamap[toupper(worldpolys$names),"color"])

