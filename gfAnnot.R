gfAnnot <- read.table("data/gf_descriptions.txt", sep = "\t", header = TRUE, as.is = TRUE)

save(list="gfAnnot", file = "data/gfAnnot.rda")
load("/Users/mikhail/Documents/Work/GenomeRunner/MDmisc/data/gfAnnot.rda")

# saveRDS(gfAnnot, file = "data/gfAnnot.Rds")
# gfAnnot <- readRDS("data/gfAnnot.Rds")
