# After running compile.popres.R, run this script in a more up-to-date
# version of R (e.g., R 3.3.2) to save the POPRES data to
# popres.RData. To reduce the file size, I provide the option of
# retaining a random subset of the SNPs on chromosomes 1-22.
# Additionally, I only retain PCs 1-4.
pcs         <- paste0("PC",1:4)
num.markers <- NULL
set.seed(1)

# Read the POPRES data from the CSV file prepared by compile.popres.R.
# This table includes the computed projection ("rotation") of the
# samples onto the PCs.
x <- read.csv("popres.csv",quote = "",comment.char = "#")
x <- transform(x,id = as.character(id))
x <- x[c("id","country","abbrv","color","longitude","latitude",pcs)]

# Load the standard deviations of the PCs.
sdev <- read.table(paste0("../original-files/data/POPRES_08_24_01.",
                          "EuroThinFinal.LD_0.8.exLD.out0-PCA.eigvals"))
sdev <- sqrt(sdev[[1]])
n    <- length(sdev)
names(sdev) <- paste0("PC",1:n)

# Load the eigenvectors ("loadings").
rotation <- read.table(paste0("../original-files/data/POPRES_08_24_01.",
                              "EuroThinFinal.LD_0.8.exLD.snpeigs"),
                       stringsAsFactors = FALSE,
                       col.names = c("marker","chr","pos",paste0("PC",1:20)))
n <- nrow(rotation)
if (is.null(num.markers)) {
  rows <- 1:n
} else {
  rows <- sort(sample(n,num.markers))
}
rotation           <- transform(rotation,chr = factor(chr))
rownames(rotation) <- rotation$marker
rotation           <- rotation[rows,c("chr","pos",pcs)]

# Save the POPRES data in an .RData file.
popres <- list(x = x,sdev = sdev[pcs],rotation = rotation)
save(list = "popres",file = "popres.RData")
