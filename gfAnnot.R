gfAnnot <- read.table("data/gf_descriptions.txt", sep = "\t", header = TRUE, as.is = TRUE)
saveRDS("data/gfAnnot.Rds")

gfAnnot <- readRDS("data/gfAnnot.Rds")
