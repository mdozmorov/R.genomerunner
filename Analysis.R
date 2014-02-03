
## ----setup, echo=F, include=FALSE, cache=F-------------------------------
library(knitr) 
opts_chunk$set(cache.path='cache/', fig.path='img/', cache=T, tidy=T, fig.keep='high', echo=F, dpi=300)
options(replace.assign=TRUE, width=65)


## ----loadLibraries, echo=TRUE, cache=FALSE-------------------------------
source("utils.R")
suppressMessages(library(Hmisc)) # For rcorr function
suppressMessages(library(gplots))
suppressMessages(library(Biobase))
suppressMessages(library(limma))


## ----loadData, echo=c(-5, -6)--------------------------------------------
# Define data subfolder to use, change to analyze different data
dname<-"data//more15_vs_tfbsEncode//"
mtx<-as.matrix(read.table(paste(dname, "matrix.txt", sep=""), sep="\t", header=T, row.names=1))
mtx<-mtx.transform(mtx) # -log10 transform p-values
# Optional: adjust columns for multiple testing. See utils.R for the function definition.
# mtx<-mtx.adjust(mtx) 


## ----preprocessData, echo=TRUE, cache=TRUE, dependson='loadData'---------
dim(mtx) # Check original dimensions
cutoff<-2 # raw p-value significance cutoff is 0.01
# What remains if we remove rows/cols with nothing significant
dim(mtx[apply(mtx, 1, function(x){sum(abs(x)>cutoff)})>0, 
        apply(mtx, 2, function(x){sum(abs(x)>cutoff)})>0])
# Trim the matrix
mtx<-mtx[apply(mtx, 1, function(x){sum(abs(x)>cutoff)})>0, 
         apply(mtx, 2, function(x){sum(abs(x)>cutoff)})>0]


## ----preprocessCorrel, echo=TRUE, dependson='preprocessData'-------------
# rcorr returns a list, [[1]] - correl coeffs, [[3]] - p-values. Type - pearson/spearman
mtx.cor<-rcorr(as.matrix(mtx), type="spearman")


## ----epigenomicVisualization, echo=c(-1, -2), fig.cap='Epigenomic similarity heatmap'----
par(oma=c(5,0,0,5)) # Adjust margins
color<-colorRampPalette(c("blue","yellow")) # Define color gradient
# Adjust clustering parameters.
# Distance: "euclidean", "maximum","manhattan", "canberra", "binary" or "minkowski".
# Clustering: "ward", "single", "complete", "average", "mcquitty", "median" or "centroid"
dist.method<-"euclidean"  
hclust.method<-"ward"
h<-heatmap.2(as.matrix(mtx.cor[[1]]), trace="none", density.info="none", col=color,distfun=function(x){dist(x, method=dist.method)}, hclustfun=function(x){hclust(x, method=hclust.method)}, cexRow=0.5, cexCol=0.5)


## ----defineClusters, echo=-1, results='hide', fig=TRUE-------------------
par(oma=c(5,0,0,5), cex=0.7)
plot(h$colDendrogram) # Plot the dendrogram only
# Cut the dentrogram into separate clusters. Tweak the height
c<-cut(h$colDendrogram,h=2.35) 
c # Check the number of clusters, and the number of members
unlink(paste(dname, "clustering.txt", sep=""))
for (i in 1:length(c$lower)){ 
  write.table(paste(i, t(labels(c$lower[[i]])), sep="\t"), paste(dname, "clustering.txt", sep=""), sep="\t", col.names=F, row.names=F, append=T)
}


## ----defineGroups, echo=TRUE, dependson='defineClusters'-----------------
eset.labels<-character() # Empty vector to hold cluster labels
eset.groups<-numeric() # Empty vector to hold cluster groups
for (i in 1:length(c$lower)) { # Go through each cluster
  # If the number of members is more than 2
  if (attr(c$lower[[i]], "members") > 2) { 
    eset.labels<-append(eset.labels, labels(c$lower[[i]]))
    eset.groups<-append(eset.groups, rep(i, length(labels(c$lower[[i]]))))
  }
}


## ----limmaOnClusters, echo=TRUE, warning=FALSE, results='hide', dependson='defineGroups'----
# Make eset out of eset.labels
eset<-new("ExpressionSet", exprs=as.matrix(mtx[, eset.labels]))
# Make model matrix
design<-model.matrix(~ 0+factor(eset.groups)) 
colnames(design)<-paste("c", unique(eset.groups), sep="")
# Create a square matrix of counts of DEGs
degs.matrix<-matrix(0, length(c$lower), length(c$lower))
# Name rows and cols by cluster name
colnames(degs.matrix)<-paste("c",seq(1,length(c$lower))); rownames(degs.matrix)<-paste("c",seq(1,length(c$lower))) 
cutoff.pval<-0.05 # p-value cutoff. Tweak
cutoff.lfc<-log(2) # log2 fold change cutoff. Tweak
for(i in colnames(design)){ 
  for(j in colnames(design)){
    # Contrasts between two clusters
    contrast.matrix<-makeContrasts(contrasts=paste(i, j, sep="-"), levels=design)
    fit <- lmFit(eset, design) 
    fit2 <- contrasts.fit(fit, contrast.matrix)
    fit2 <- eBayes(fit2)
    degs<-topTable(fit2, number=dim(exprs(eset))[[1]], adjust.method="BH", p.value=cutoff.pval, lfc=cutoff.lfc)
    print(paste(i, "vs.", j, ", number of degs:", dim(degs)[[1]]))
    # Keep the number of DEGs in the matrix
    degs.matrix[as.numeric(sub("c","",i)), as.numeric(sub("c","",j))]<-dim(degs)[[1]]
    if(dim(degs)[[1]]>0) {
      # Average values in clusters i and j
      i.av<-rowMeans(matrix(exprs(eset)[rownames(degs),eset.groups == as.numeric(sub("c","",i))], nrow=dim(degs)[[1]]))
      j.av<-rowMeans(matrix(exprs(eset)[rownames(degs),eset.groups == as.numeric(sub("c","",j))], nrow=dim(degs)[[1]]))
      i.vs.j<-rep(paste(i,"vs.",j), dim(degs)[[1]])
      # Put it all together in a file
      write.table(cbind(degs, i.vs.j, i.av , j.av), paste(dname, "degs.txt", sep=""), sep="\t", col.names=NA, append=T)
    }
  }
}


## ----maxminCorr, echo=TRUE, eval=TRUE, results='hide', dependson='preprocessCorrel'----
mtx.cor1<-mtx.cor[[1]]
# We don't need to consider perfect correlations, zero them out
diag(mtx.cor1)<-0
# Print top correlated parameters on screen
for (i in head(unique(mtx.cor1[order(mtx.cor1,decreasing=T)]))) {print(which(mtx.cor1 == i, arr.ind=T))}
unlink(paste(dname, "maxmin_correlations.txt", sep=""))
for (i in 1:ncol(mtx.cor1)) write.table(paste(colnames(mtx.cor1)[i],"correlates with",
                                              colnames(mtx.cor1)[which(mtx.cor1[i,] == max(mtx.cor1[i,]))], 
                                              "at corr. coeff.",formatC(mtx.cor1[i,which(mtx.cor1[i,] == max(mtx.cor1[i,]))]),
                                              "anticorrelates with",
                                              colnames(mtx.cor1)[which(mtx.cor1[i,] == min(mtx.cor1[i,]))],
                                              "at corr. coeff.",formatC(mtx.cor1[i,which(mtx.cor1[i,] == min(mtx.cor1[i,]))]),sep="|"),
                                        paste(dname, "maxmin_correlations.txt", sep=""), append=T, sep="\t", col.names=F, row.names=F) 


## ----enrichmentCutoffs, echo=TRUE, eval=TRUE, results='hide', fig.show='hold', dependson='preprocessCorrel'----
 # Define minimum number of times a row/col should have values above the cutoffs
numofsig<-1
dim(mtx) # Original dimensions
# Check summary and set p-value and variability cutoffs as 3rd quantiles of their distributions
summary(as.vector(abs(mtx))); cutoff.p<-summary(as.vector(abs(mtx)))[[5]]
summary(as.vector(apply(abs(mtx),1,sd))); cutoff.sd<-summary(as.vector(apply(abs(mtx),1,sd)))[[5]]
# Check visual distributions and set p-value and variability cutoffs manually
hist(as.vector(mtx), breaks=20, main="Distribution of -log10-transformed p-values", xlab="-log10-transformed p-values")
hist(c(as.vector(apply(mtx,1,sd)), as.vector(apply(mtx,2,sd))), breaks=20, main="Distribution of SD across rows and columns", xlab="SD")
cutoff.p<-4; cutoff.sd<-0.5


## ----trimCutoffs, echo=TRUE, results='hide', dependson='enrichmentCutoffs'----
mtx.gf<-mtx
# Remove rows/cols that do not show significant p-values and variability less than numofsig times
mtx.gf<-mtx.gf[apply(mtx.gf, 1, function(row){sum(abs(row)>cutoff.p)>=numofsig}), 
               apply(mtx.gf, 2, function(col){sum(abs(col)>cutoff)>=numofsig})]
mtx.gf<-mtx.gf[apply(mtx.gf, 1, sd)>cutoff.sd,
               apply(mtx.gf, 2, sd)>cutoff.sd]
dim(mtx.gf) # Dimensions after trimming


## ----enrichmentVisualization, echo=c(-8, -9, -10, -11, -12, -14), results='hide', warning=FALSE, fig=TRUE, fig.cap='Heatmap of the enrichment results', fig.show='asis', dependson='trimCutoffs', dependson='enrichmentCutoffs'----
# Adjust clustering parameters.
# Distance: "euclidean", "maximum","manhattan", "canberra", "binary" or "minkowski".
# Clustering: "ward", "single", "complete", "average", "mcquitty", "median" or "centroid"
dist.method<-"maximum"
hclust.method<-"ward"
granularity=7
# my.breaks<-c(seq(min(mtx.gf), -cutoff.p, length.out=granularity), seq(cutoff.p, max(mtx.gf), length.out=granularity)) 
# my.breaks<-c(seq(min(mtx.gf), -2, length.out=granularity), seq(2, max(mtx.gf), length.out=granularity))
my.breaks<-c(seq(min(mtx.gf), max(mtx.gf)))
par(oma=c(10,0,0,3))
color<-colorRampPalette(c("blue", "yellow"))
h<-heatmap.2(as.matrix(mtx.gf), distfun=function(x){dist(x,method=dist.method)}, hclustfun=function(x){hclust(x,method=hclust.method)}, dendrogram="none", breaks=my.breaks,  col=color, lwid=c(1.5,3), lhei=c(1.5,4), key=T,  symkey=T, keysize=0.01, density.info="none", trace="none",  cexCol=1.0, cexRow=0.8)


