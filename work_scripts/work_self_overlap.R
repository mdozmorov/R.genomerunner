mtx<-as.matrix(read.table("overlap.txt",sep="\t",header=T,row.names=1))
mtx[mtx == 1] <- 0.9999
diag(mtx)<-1
dist.method<-"maximum"  # "euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski"
hclust.method<-"ward" # "ward", "single", "complete", "average", "mcquitty", "median" or "centroid"
color<-colorRampPalette(c("blue", "red"))
color<-colorRampPalette(c("blue", "yellow"))
# To have p-values overlayed on the cells: cellnote=formatC(mtx, digits=2)
heatmap.2(as.matrix(mtx), distfun=function(x){dist(x,method=dist.method)}, hclustfun=function(x){hclust(x,method=hclust.method)}, dendrogram="column", col=color, lwid=c(1.5,3), lhei=c(1.5,4), key=T,  symkey=T, keysize=0.01, density.info="none", trace="none",  cexCol=1.0, cexRow=0.8)#, cellnote=formatC(mtx, digits=2), 

heatmap.2(as.matrix(mtx), Rowv=F, Colv=F, dendrogram="column", col=color, lwid=c(1.5,3), lhei=c(1.5,4), key=T,  symkey=F, keysize=0.01, density.info="none", trace="none",  cexCol=1.0, cexRow=0.8, cellnote=formatC(mtx, digits=2))