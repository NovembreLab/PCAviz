# This short script demonstrates plotting of samples projected onto
# principal components, after applying a rotation to the projection,
# in which samples are coloured based on an expert-provided population
# assignment. Additionally, country-level summary statistics (medians)
# are visualized in the same PC plot, and these statistics are
# visually compared against a map of Europe. This code closely
# reproduces Fig. 1a in the "Genes mirror geography" paper.
#
# This script has been tested in R version 2.5.1, with sp package
# version 0.9-28. Note: since the map rendering can be very slow, I
# recommend using ThinLinc instead of X11 forwarding to display
# graphics.
#
library(sp)
source("../functions.R")

# SCRIPT PARAMETERS
# -----------------
x  <- 1     # First PC ("PC1").
y  <- 2     # Second PC ("PC2").
sx <- 1     # Scaling factor for PC1.
sy <- (-1)  # Scaling factor for PC2.

# Load the PCA results.
PCA <- read.table(paste("../../original-files/data/POPRES_08_24_01.",
                        "EuroThinFinal.LD_0.8.exLD.out0-PCA.eigs",sep=""),
                  col.names = c("FamID","ID",paste("PC",1:20,sep=""),"flag"))

# Load the POPRES population and geographic annotations.
ddata <- read.table("../../original-files/data/PCA.txt",sep = "\t",
                    header = TRUE,stringsAsFactors = FALSE)

# Fix up some of the country-of-origin annotations for the plots.
poslabel    <- fix.country.names(ddata$plabels)
collabelmap <- poslabel

# Load the map-color annotations.
coldatamap <- read.table("../../original-files/data/ColorTablePCmap.txt",
                         sep = "\t",stringsAsFactors = FALSE,
                         col.names = c("cntry","color"))
row.names(coldatamap) <- toupper(coldatamap$cntry)

# Get the abbreviated population labels assigned to each of the POPRES
# samples.
abbrvlabel <- abbrvlabels(ddata$plabels)

# Retrieve the geographic co-ordinate data for the countries.
load("../../original-files/data/world_countries.rda")
xvals <- coordinates(world_countries)[match(toupper(levels(factor(poslabel))),
                                            world_countries$names),1]
yvals <- coordinates(world_countries)[match(toupper(levels(factor(poslabel))),
                                            world_countries$names),2]
posdata            <- data.frame(Long = xvals,Lat = yvals)
row.names(posdata) <- levels(factor(poslabel))

# Get the labels, colors and geographic co-ordinates assigned to each
# of the POPRES samples.
rows    <- match(PCA$ID,ddata$ID)
alabels <- abbrvlabel[rows]
clabels <- collabelmap[rows]
plabels <- poslabel[rows]
col     <- coldatamap[toupper(clabels),"color"]
xvals   <- posdata[as.character(plabels),"Long"]
yvals   <- posdata[as.character(plabels),"Lat"]

# Determine the rotation angle (phi) that offers the best fit between the
# PCs and the geographic co-ordinates.
phivals <- seq(0,2*pi,length = 360)
calccor <- function (p)
  cor(rotPC1(p,sx,sy),yvals) + cor(rotPC2(p,sx,sy),xvals)
corvals <- apply(as.matrix(phivals),1,calccor)
phi <- phivals[which.max(corvals)]
ang <- round(phi/(2*pi)*360)
cat("phi =",phi,"\n")

# Create a new graphics window (this command may not work on all devices).
x11()

plot(0,0,xlab = "",ylab = "",type = "n",asp = 1,frame.plot = FALSE,
     xaxt = "n",yaxt = "n",xlim = c(-0.07,0.09))

# Some plotting parameters.
x0    <- 0
y0    <- 0
m     <- 0.05
caxis <- "black"

# Draw the arrows.
m <- 0.095
arrows(0,0,x0 + m*cos(-phi),y0 + m*sin(-phi),col = caxis,length = 0.1)
m <- 0.07
arrows(0,0,x0 + m*cos(-phi + pi/2),y0 + m*sin(-phi + pi/2),col = caxis,
       length = 0.1)
m <- 0.06
arrows(0,0,x0 - m*cos(-phi),y0 - m*sin(-phi),col = caxis,length = 0.1)
m <- 0.0725
arrows(0,0,x0 - m*cos(-phi + pi/2),y0 - m*sin(-phi + pi/2),col = caxis,
       length = 0.1)

# Plot the POPRES samples based on their projection onto PCs 1 and 2,
# after the specified rotation; each sample is represented by the
# abbreviated country-of-origin code.
i         <- sample(length(PCA$ID))
text(rotPC2(phi,sx,sy)[i],rotPC1(phi,sx,sy)[i],alabels[i],col = col[i],
     cex = 0.75)

# Next, plot the median lat-long locations by country, with some
# manual adjustment of the final locations.
PC2median <- tapply(rotPC2(phi,sx,sy),alabels,median)
PC1median <- tapply(rotPC1(phi,sx,sy),alabels,median)
colmedian <- tapply(col,alabels,function(x) x[1])

PC1medianOrig <- PC1median
PC2medianOrig <- PC2median
PC2median[names(PC1median)=="BG"] <- PC2median[names(PC1median)=="BG"]-0.0055
PC2median[names(PC1median)=="MK"] <- PC2median[names(PC1median)=="MK"]+0.005
PC1median[names(PC1median)=="SI"] <- PC1median[names(PC1median)=="SI"]-0.005
PC2median[names(PC1median)=="HR"] <- PC2median[names(PC1median)=="HR"]+0.005
PC1median[names(PC1median)=="RS"] <- PC1median[names(PC1median)=="RS"]-0.005
PC2median[names(PC1median)=="BA"] <- PC2median[names(PC1median)=="BA"]+0.005
PC1median[names(PC1median)=="CY"] <- PC1median[names(PC1median)=="CY"]-0.0055

segments(PC2median,PC1median,PC2medianOrig,PC1medianOrig,col = 1)
points(PC2medianOrig,PC1medianOrig,col = colmedian,pch = 16,cex = 1.5)
points(PC2median,PC1median,col = colmedian,pch = 16,cex = 4) 
text(PC2median,PC1median,names(PC1median))

# Create a new graphics window (this command may not work on all devices).
x11()

# Draw a map of Europe with countries coloured in the same way as they
# were in the first plot.
plot(world_countries,xlim = c(-8,35),ylim = c(35,60),
     col = coldatamap[toupper(world_countries$names),"color"])

# Get the country abbreviations.
r         <- tapply(plabels,factor(alabels),function(x) x[1])
colmedian <- tapply(col,alabels,function(x) x[1])

# Make a few small adjustments to the co-ordinates for plotting.
xleg <- posdata[r,1]
yleg <- posdata[r,2]
xleg[r=="Russian Federation"]=xleg[r=="Russian Federation"]-2
xleg[r=="Albania"]= xleg[r=="Albania"]-1
xleg[r=="Slovenia"]=xleg[r=="Slovenia"]-1
yleg[r=="Slovenia"]=yleg[r=="Slovenia"]-1
yleg[r=="Kosovo"]=yleg[r=="Kosovo"]+0.25
xleg[r=="Kosovo"]=xleg[r=="Kosovo"]+0.35
xleg[r=="Bosnia and Herzegovina"]=xleg[r=="Bosnia and Herzegovina"]-1
xleg[r=="Yugoslavia"]=xleg[r=="Serbia and Montenegro"]-1.25
yleg[r=="Yugoslavia"]=yleg[r=="Serbia and Montenegro"]
yleg[r=="Serbia and Montenegro"]=yleg[r=="Serbia and Montenegro"]+1.25

# Draw the abbreviated country labels.
points(xleg,yleg,col = colmedian,pch = 16,cex = 2.95)
text(xleg,yleg,levels(factor(alabels)),cex = 0.75)

