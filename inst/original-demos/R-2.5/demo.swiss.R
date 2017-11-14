# This short script demonstrates plotting of samples projected onto
# principal components, after applying a rotation to the projection,
# in which samples are coloured based on an expert-provided population
# assignment. This code roughly reproduces Fig. 1b in the "Genes
# mirror geography" paper.
#
# This script has been tested in R version 2.5.1.
#
source("../functions.R")

# SCRIPT PARAMETERS
# -----------------
sx  <- 1     # Scaling factor for PC1.    
sy  <- (-1)  # Scaling factor for PC2.
phi <- -16   # Angle of rotation, in radians (equal to 160 degrees).

# Load the PCA results.
PCA <- read.table(paste("../../original-files/data/POPRES_08_24_01.",
                        "EuroThinFinal.LD_0.8.exLD.out0-PCA.eigs",sep=""),
                  col.names = c("FamID","ID",paste("PC",1:20,sep=""),"flag"))

# Load the POPRES population and geographic annotations.
ddata <- read.table("../../original-files/data/PCA.txt",sep = "\t",
                    header = TRUE,stringsAsFactors = FALSE)

# Get the abbreviated population labels assigned to each of the POPRES
# samples.
rows    <- match(PCA$ID,ddata$ID)
alabels <- abbrvlabels.swiss(ddata$plabels[rows])

# Plot the POPRES samples projected onto PCs 1 and 2, after applying a
# rotation, coloured by their population label.
plot(rotPC2(phi,sx,sy)[!is.element(alabels,c("CH-F","CH-I","CH-G"))],
     rotPC1(phi,sx,sy)[!is.element(alabels,c("CH-F","CH-I","CH-G"))],
     xlab = "",ylab = "",asp = 1,pch = 16,cex = 0.5,col = "lightgrey",
     xlim = c(-0.03,0.03),ylim = c(-0.03,0.03))
points(rotPC2(phi,sx,sy)[alabels == "IT"],
       rotPC1(phi,sx,sy)[alabels == "IT"],
       col = "lightblue",cex = 1,pch = 17)
points(rotPC2(phi,sx,sy)[alabels == "DE"],
       rotPC1(phi,sx,sy)[alabels=="DE"],
       col="pink",cex = 1.2,pch = 18)
points(rotPC2(phi,sx,sy)[alabels == "FR"],
       rotPC1(phi,sx,sy)[alabels == "FR"],
       col = "greenyellow",cex = 1,pch = 15)
points(rotPC2(phi,sx,sy)[alabels == "CH-I"],
       rotPC1(phi,sx,sy)[alabels == "CH-I"],
       col = "blue",cex = 1.25,pch = 2)
points(rotPC2(phi,sx,sy)[alabels == "CH-G"],
       rotPC1(phi,sx,sy)[alabels == "CH-G"],
       col = "red",cex = 1.25,pch = 5)
points(rotPC2(phi,sx,sy)[alabels == "CH-F"],
       rotPC1(phi,sx,sy)[alabels == "CH-F"],
       col = "green4",cex = 1.25,pch = 22)

