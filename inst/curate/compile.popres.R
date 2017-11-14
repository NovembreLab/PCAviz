# Compile the POPRES data into a convenient form. The data preparation
# steps implemented here are drawn from demo.biplot.R. This script has
# been tested in R version 2.5.1, with sp package version 0.9-28.
library(sp)
source("../demo/functions.R")

# SCRIPT PARAMETERS
# -----------------
x  <- 1     # First PC ("PC1").
y  <- 2     # Second PC ("PC2").
sx <- 1     # Scaling factor for PC1.
sy <- (-1)  # Scaling factor for PC2.

# Load the PCA results.
PCA <- read.table(paste("../original-files/data/POPRES_08_24_01.",
                        "EuroThinFinal.LD_0.8.exLD.out0-PCA.eigs",sep=""),
                  col.names = c("FamID","ID",paste("PC",1:20,sep=""),"flag"))

# Load the POPRES population and geographic annotations.
ddata <- read.table("../original-files/data/PCA.txt",sep = "\t",
                    header = TRUE,stringsAsFactors = FALSE)

# Fix up some of the country-of-origin annotations for the plots.
poslabel    <- fix.country.names(ddata$plabels)
collabelmap <- poslabel

# Load the map-color annotations.
coldatamap <- read.table("../original-files/data/ColorTablePCmap.txt",
                         sep = "\t",stringsAsFactors = FALSE,
                         col.names = c("cntry","color"))
row.names(coldatamap) <- toupper(coldatamap$cntry)

# Retrieve the geographic co-ordinate data for the countries.
load("../original-files/data/world_countries.rda")
xvals <- coordinates(world_countries)[match(toupper(levels(factor(poslabel))),
                                            world_countries$names),1]
yvals <- coordinates(world_countries)[match(toupper(levels(factor(poslabel))),
                                            world_countries$names),2]
posdata            <- data.frame(Long = xvals,Lat = yvals)
row.names(posdata) <- levels(factor(poslabel))

# Get the abbreviated population labels assigned to each of the POPRES
# samples.
abbrvlabel <- abbrvlabels(ddata$plabels)

# Get the labels, colors and geographic co-ordinates assigned to each
# of the POPRES samples.
rows    <- match(PCA$ID,ddata$ID)
alabels <- abbrvlabel[rows]
plabels <- poslabel[rows]
col     <- coldatamap[toupper(plabels),"color"]
xvals   <- posdata[as.character(plabels),"Long"]
yvals   <- posdata[as.character(plabels),"Lat"]

# Fix the Macedonia annotations so that they don't contain commas.
plabels[grep("Macedonia",plabels)] <- "Macedonia"

# Combine the data into one large table.
combined.data <-
  cbind(PCA["ID"],
        data.frame(country   = plabels,
                   abbrv     = alabels,
                   color     = col,
                   longitude = round(xvals,digits = 4),
                   latitude  = round(yvals,digits = 4)),
        PCA[paste("PC",1:20,sep="")])

# Write the POPRES data to a comma-delimited file.
write.csv(combined.data,"popres.csv",quote = FALSE,row.names = FALSE)
