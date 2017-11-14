# This script creates scatterplots comparing the median projection, by
# country-of-origin, of the POPRES samples onto principal components
# (PC) 1 and 2, against geographic co-ordinates (latitude and
# longitude). Compare these plots against Figure 1 (b) and (c) in the
# supplementary materials of the "Genes mirror geography" paper.
#
# This script has been tested in R version 2.5.1, with sp package
# version 0.9-28.
#
library(sp)
source("../functions.R")

# SCRIPT PARAMETERS
# -----------------
sx    <- (-1)   # Scaling factor for PC1.
sy    <- 1      # Scaling factor for PC2.
phi   <- (-16)  # Angle (in degrees).
alpha <- 0.025  # Defines upper and lower quantiles in plots.

# Convert angle from degrees to radians.
phi <- 2*pi*(phi/360)

# Load the PCA results.
PCA <- read.table(paste("../../original-files/data/POPRES_08_24_01.",
                        "EuroThinFinal.LD_0.8.exLD.out0-PCA.eigs",sep=""),
                  col.names = c("FamID","ID",paste("PC",1:20,sep=""),"flag"))

# Load the POPRES population and geographic annotations.
ddata <- read.table("../../original-files/data/PCA.txt",sep = "\t",
                    header = TRUE,stringsAsFactors = FALSE)

# Fix up some of the country-of-origin annotations for drawing the maps.
poslabel    <- fix.country.names(ddata$plabels)
maplabel    <- fix.country.names.maps(ddata$plabels)
collabelmap <- poslabel

# Retrieve the geographic co-ordinate data for the countries.
load("../../original-files/data/world_countries.rda")
xvals <- coordinates(world_countries)[match(toupper(levels(factor(poslabel))),
                                            world_countries$names),1]
yvals <- coordinates(world_countries)[match(toupper(levels(factor(poslabel))),
                                            world_countries$names),2]
posdata            <- data.frame(Long = xvals,Lat = yvals)
row.names(posdata) <- levels(factor(poslabel))

# Load the map-color annotations.
coldatamap <- read.table("../../original-files/data/ColorTablePCmap.txt",
                         sep = "\t",stringsAsFactors = FALSE,
                         col.names = c("cntry","color"))
row.names(coldatamap) <- toupper(coldatamap$cntry)

# Get the abbreviated population labels assigned to each of the POPRES
# samples.
abbrvlabel <- abbrvlabels(ddata$plabels)

# Get the labels and colors assigned to the POPRES samples.
rows     <- match(PCA$ID,ddata$ID)
maplabel <- maplabel[rows]
alabels  <- abbrvlabel[rows]
plabels  <- poslabel[rows]
clabels  <- collabelmap[rows]

# Compute the median and standard error (s.e.) for each country label.
se <- function (x)
  sd(x)/sqrt(length(x))
PC2median <- tapply(rotPC2(phi,sx,sy),maplabel,median)
PC1median <- tapply(rotPC1(phi,sx,sy),maplabel,median)
PC2se     <- tapply(rotPC2(phi,sx,sy),maplabel,se)
PC1se     <- tapply(rotPC1(phi,sx,sy),maplabel,se)
labels    <- tapply(plabels,factor(maplabel),function(x) x[1])

# Get the latitude and longitude for each of the country labels.
i    <- tapply(plabels,factor(maplabel),function (x) x[1])
xleg <- posdata[i,1]
yleg <- posdata[i,2]    

# Get the color and abbreviated label assigned to each country. 
col       <- coldatamap[toupper(clabels),"color"]
colmedian <- tapply(col,maplabel,function(x) x[1])
labelsA   <- tapply(alabels,maplabel,function(x) x[1])

# Create a scatterplot comparing PC1 (after a rotation) against latitude.
plot(yleg,PC1median,pch = 16,type = "n",xlab = "Latitude",
     ylab = "South-North in rotated PC1-PC2 space")
points(yleg,PC1median,pch = 16,cex = 3,col = colmedian)
text(yleg,PC1median,labelsA,cex = 0.75)

# Plot lines representing the upper and lower tails of the residuals
# under a simple linear regression.
lPC1      <- lm(PC1median ~ yleg)
residsPC1 <- abs(lPC1$residuals)
abline(lPC1)
abline(lPC1$coefficients[1] + quantile(residsPC1,1-alpha),
       lPC1$coefficients[2],lty = 2)
abline(lPC1$coefficients[1] - quantile(residsPC1,1 - alpha),
       lPC1$coefficients[2],lty = 2)

# Create a scatterplot comparing PC2 (after a rotation) against longitude.
x11()
plot(xleg,PC2median,pch = 16,type = "n",xlab = "Longitude",
     ylab = "West-East in rotated PC1-PC2 space")
points(xleg,PC2median,pch = 16,cex = 3,col = colmedian)
text(xleg,PC2median,labelsA,cex = 0.75)

# Plot lines representing the upper and lower tails of the residuals
# under a simple linear regression.
lPC2      <- lm(PC2median ~ xleg)
residsPC2 <- abs(lPC2$residuals)
abline(lPC2)
abline(lPC2$coefficients[1] + quantile(residsPC2,1 - alpha),
       lPC2$coefficients[2],lty = 2)
abline(lPC2$coefficients[1] - quantile(residsPC2,1 - alpha),
       lPC2$coefficients[2],lty = 2)
