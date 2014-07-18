source("MINE.r")
# === convert matrix to MINE format
mtx <- as.data.frame(read.table("clipboard", sep="\t", header=T, row.names=1)) #Matrix of unadjusted p-values
write.table(t(mtx),"gwascatalog_more50_vs_genes_fisher.csv",sep=",",row.names=T,col.names=F) #Transpose, do not output header

MINE("gwascatalog_more50_vs_genes_fisher.csv","all.pairs")

