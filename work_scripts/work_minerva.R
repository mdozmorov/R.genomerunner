install.packages("minerva")
library(minerva)
library(gplots)
library(RColorBrewer)
display.brewer.all(n=NULL, type="all", select=NULL, exact.n=TRUE)
par(oma=c(5,0,0,5)) #Make right and bottom margins larger

mtx<-read.table("gwascatalog_more50_vs_genes_fisher.txt",sep="\t",row.names=1, header=T)
mtx<-mtx[apply(mtx,1,sd)>=0.5,] # Remove rows with SD=0
mtx<-mtx[,apply(mtx,2,sd)>=0.5] # Remove cols with SD=0
IAC<-mine(mtx)$MIC
diag(IAC) # Stopped when diag not always = 1
names(m)
dim(m$MIC)
IAC<-m$MIC
colnames(IAC)<-colnames(mtx)
rownames(IAC)<-colnames(mtx)
hist(IAC,n=50)
median(IAC)
# Quantile normalization of the whole matrix
C<-IAC
C[order(C)] <- 1:length(C)/length(C) 
IAC<-C

granularity=6
dist.method<-"euclidean"  # "euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski"
hclust.method<-"ward" # "ward", "single", "complete", "average", "mcquitty", "median" or "centroid"
my.breaks<-c(seq(min(IAC),median(IAC),length.out=granularity),seq(median(IAC),max(IAC),length.out=granularity))
# With breaks, as is data
heatmap.2(as.matrix(IAC), distfun=function(x){dist(x,method=dist.method)}, hclustfun=function(x){hclust(x,method=hclust.method)}, dendrogram="none", breaks=my.breaks, col=rev(brewer.pal(2*granularity-1,'RdBu')),lwid=c(1.5,3), lhei=c(1.5,4), key=T,  symkey=F, keysize=0.01, density.info="none", trace="none",  cexCol=1.0, cexRow=1)#, cellnote=formatC(mtx, digits=2), notecex=1, notecol='yellow') 

heatmap.2(t(scale(t(scale(mtx)))), distfun=function(x){dist(x,method=dist.method)}, hclustfun=function(x){hclust(x,method=hclust.method)}, dendrogram="none", col=rev(brewer.pal(11,'RdBu')),lwid=c(1.5,3), lhei=c(1.5,4), key=T,  symkey=T, keysize=0.01, density.info="nomne", trace="none",  cexCol=1.0, cexRow=1)#, cellnote=formatC(mtx, digits=2), notecex=1, notecol='yellow') 


# Standardizind the data
mtx<-scale(mtx,center=T,scale=T) # Standardize columns
