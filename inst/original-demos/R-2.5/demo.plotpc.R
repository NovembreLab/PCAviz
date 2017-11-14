# This script demonstrates colouring of countries in a map according
# to their projection onto a principal component (PC).
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
x  <- 1  # Show results for this principal component.
sx <- 1  # Scaling factor for selected PC.

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

# Retrieve the geographic co-ordinate data for the countries.
load("../../original-files/data/world_countries.rda")
xvals <- coordinates(world_countries)[match(toupper(levels(factor(poslabel))),
                                            world_countries$names),1]
yvals <- coordinates(world_countries)[match(toupper(levels(factor(poslabel))),
                                            world_countries$names),2]

# Get POPRES sample annotations.
rows    <- match(PCA$ID,ddata$ID)
clabels <- collabelmap[rows]

# Retrieve all country labels with >1 sample.
aboveThresh <- levels(factor(clabels))[table(clabels) >= 2]
labelsThin  <- clabels[is.element(clabels,aboveThresh)]

# Compute mean and standard deviation of the PC projection by country label.
PCvals   <- sx*PCA[is.element(clabels,aboveThresh),paste("PC",x,sep="")]
meanvals <- tapply(PCvals,factor(labelsThin),median)
sdvals   <- tapply(PCvals,factor(labelsThin),sd)

# Color each country according to its mean PC value: smallest (most
# negative) values are colored red; largest (most positive) values are
# coloured pale yellow. Countries with insufficient data (NA) are
# shown in light gray ("gainsboro").
bins   <- seq(1.01 * min(PCvals),1.01 * max(PCvals),length = 15)
colors <- rainbow(length(bins) - 1,s = 0.6,start = 0/6,end = 4/6)
colors <- heat.colors(length(bins) - 1)
countrycolors <- apply(meanvals,1,function (x) {
  for(i in 1:length(bins))
    if (x < bins[i]) {
      return(i-1)
      break
    }
  return(0)
})
y <- colors[countrycolors[match(world_countries$names,
                                toupper(names(countrycolors)))]]
y[is.na(y)] <- "gainsboro"

# Plot the map with the countries colored according to their
# projection onto the selected PC.
plot(world_countries,xlim = c(-12,55),ylim = c(35,65),col = y,
     main = paste("PC",x,sep=""))
