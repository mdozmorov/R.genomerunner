library(dplyr)
# Read in gf_descriptions file created by gr-dbcreator 
gf_descriptions_nodups_hg19 <- read.table("data/db_6.00_02-07-2016/grsnp_db/hg19/gf_descriptions_nodups_hg19.txt.gz", sep = "\t", header = T, stringsAsFactors = F)
# Sort it
gf_descriptions_nodups_hg19 <- gf_descriptions_nodups_hg19[order(gf_descriptions_nodups_hg19[, "category"],
                                                                 gf_descriptions_nodups_hg19[, "cell"],
                                                                 gf_descriptions_nodups_hg19[, "factor"]), ]
# Available ENCODE cells
gf_ENCODE_cells <- gf_descriptions_nodups_hg19$cell[ grep("^E", gf_descriptions_nodups_hg19$cell, invert = T)] %>% unique
# Read in ENCODE cell types
ENCODE_cells <- read.table("data/ENCODE_cells.txt", sep = "\t", header = T, quote = "\"", fill = T, stringsAsFactors = F)
tail(sort(table(ENCODE_cells$cell)))
setdiff(toupper(gf_ENCODE_cells), toupper(ENCODE_cells$cell))


# Convert descriptions to RDA object
gfAnnot <- read.table("data/gf_descriptions.txt", sep = "\t", header = TRUE, as.is = TRUE)

save(list="gfAnnot", file = "data/gfAnnot.rda")
load("/Users/mikhail/Documents/Work/GenomeRunner/MDmisc/data/gfAnnot.rda")

# saveRDS(gfAnnot, file = "data/gfAnnot.Rds")
# gfAnnot <- readRDS("data/gfAnnot.Rds")
