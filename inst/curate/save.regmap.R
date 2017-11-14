# After running compile.regmap.R, run this script in a more up-to-date
# version of R (e.g., R 3.3.2) to save the RegMap data to regmap.RData.
# We retain the first 10 PCs.
pcs <- paste0("PC",1:10)

# Read the RegMap data from the CSV file prepared by compile.regmap.R.
regmap <- read.csv("regmap.csv",quote = "",comment.char = "#",
                   as.is = c("nativename","firstname","surname","site",
                             "region"))
regmap <- transform(regmap,
                    array_id = as.character(array_id),
                    ecotype_id = as.character(ecotype_id))                    
cols   <- c(names(regmap)[1:11],pcs)
regmap <- regmap[cols]

# If the region is the empty string, set it to missing, then convert
# the "region" column to a factor.
rows                  <- which(regmap$region == "")
regmap[rows,"region"] <- NA
regmap$region         <- factor(regmap$region)

# Save the RegMap data in an .RData file.
save(list = "regmap",file = "regmap.RData")
