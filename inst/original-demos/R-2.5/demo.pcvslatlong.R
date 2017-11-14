# This script creates scatterplots comparing the projection of the
# POPRES samples onto principal components (PC) 1 and 2 against
# geographic co-ordinates (latitude and longitude).
#
# This script has been tested in R version 2.5.1, with sp package
# version 0.9-28.
#
library(sp)
source("../functions.R")

# SCRIPT PARAMETERS
# -----------------
sx  <- (-1)   # Scaling factor for PC1.
sy  <- (-1)   # Scaling factor for PC2.
phi <- (-16)  # Angle (in degrees).

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

# Get the labels and colors assigned to the POPRES samples.
rows    <- match(PCA$ID,ddata$ID)
alabels <- abbrvlabel[rows]
plabels <- poslabel[rows]
clabels <- collabelmap[rows]

# Retrieve the geographic co-ordinate data for the countries.
load("../../original-files/data/world_countries.rda")
xvals <- coordinates(world_countries)[match(toupper(levels(factor(poslabel))),
                                            world_countries$names),1]
yvals <- coordinates(world_countries)[match(toupper(levels(factor(poslabel))),
                                            world_countries$names),2]
posdata            <- data.frame(Long = xvals,Lat = yvals)
row.names(posdata) <- levels(factor(poslabel))

# Get the geographic co-ordinates and colors for the POPRES samples.
xvals <- posdata[as.character(plabels),"Long"]
yvals <- posdata[as.character(plabels),"Lat"]
col   <- coldatamap[toupper(clabels),"color"]

# Convert angle from degrees to radians.
phi <- 2*pi*(phi/360)

# Estimate the linear relationship between PC1 (after rotation) and
# latitude.
lPC1      <- lm(rotPC1(phi,sx,sy) ~ yvals)
residsPC1 <- abs(lPC1$residuals)
alpha     <- 0.025
r         <- abs(resid(lPC1))
outlier   <- r > quantile(r,0.975)

# Create a scatterplot to visualize the relationship between PC1 and
# latitude.
plot(yvals,rotPC1(phi,sx,sy),xlab = "Latitude",
     ylab = "North-south in rotated PC1-PC2 space",pch = 16,type = "n")
points(yvals[!outlier],rotPC1(phi,sx,sy)[!outlier],pch = 16,
       col = col[!outlier])
text(yvals[outlier],rotPC1(phi,sx,sy)[outlier],alabels[outlier],
     cex = 0.9,col = col[outlier])

abline(lPC1)
abline(lPC1$coefficients[1] + quantile(residsPC1,1 - alpha),
       lPC1$coefficients[2],lty = 2)
abline(lPC1$coefficients[1] - quantile(residsPC1,1 - alpha),
       lPC1$coefficients[2],lty = 2)

# Calculate the correlation between PC1 and latitude, with and without
# outliers.
print(paste("Lat Correlation = ",signif(cor(yvals,rotPC1(phi,sx,sy)),3)))
print(paste("Lat Correlation (no outliers) = ",
            signif(cor(yvals[!outlier],rotPC1(phi,sx,sy)[!outlier]),3)))

# Estimate the linear relationship between PC1 (after rotation) and
# latitude.
lPC2      <- lm(rotPC2(phi,sx,sy) ~ xvals)
residsPC2 <- abs(lPC2$residuals)
alpha     <- 0.025
r         <- abs(resid(lPC2))
outlier   <- r > quantile(r,0.975)

# Create a scatterplot to visualize the relationship between PC2 and
# longitude.
x11()
plot(xvals,rotPC2(phi,sx,sy),xlab = "Longitude",
     ylab = "East-west in rotated PC1-PC2 space",type = "n")
points(xvals[!outlier],rotPC2(phi,sx,sy)[!outlier],pch = 16,
       col = col[!outlier])
text(xvals[outlier],rotPC2(phi,sx,sy)[outlier],alabels[outlier],
     cex = 0.9,col = col[outlier])

abline(lPC2)
abline(lPC2$coefficients[1] + quantile(residsPC2,1 - alpha),
       lPC2$coefficients[2],lty = 2)
abline(lPC2$coefficients[1] - quantile(residsPC2,1 - alpha),
       lPC2$coefficients[2],lty = 2)

# Calculate the correlation between PC2 and longitude, with and without
# outliers.
print(paste("Long Correlation = ",signif(cor(xvals,rotPC2(phi,sx,sy)),3)))
print(paste("Long Correlation (no outliers) = ",
            signif(cor(xvals[!outlier],rotPC2(phi,sx,sy)[!outlier]),3)))
