
## ----setup, echo=F, include=FALSE, cache=F-------------------------------
library(knitr) 
opts_chunk$set(cache.path='cache/', fig.path='img/', cache=T, tidy=T, fig.keep='high', echo=F, dpi=300)
options(replace.assign=TRUE, width=65)
#  listing <- function(x, options) {
#    paste("\\begin{lstlisting}[basicstyle=\\ttfamily,breaklines=true]\n",
#      x, "\\end{lstlisting}\n", sep = "")
#  }
#  knit_hooks$set(source=listing, output=listing)


## ----loadLibraries, echo=TRUE, cache=FALSE-------------------------------
source("utils.R")
suppressMessages(library(Hmisc)) # For rcorr function
suppressMessages(library(gplots))
suppressMessages(library(Biobase))
suppressMessages(library(limma))


## ----loadData, echo=c(-7, -8, -10, -11)----------------------------------
# Define output and data subfolders to use, change to analyze different data
rname<-"data//results//" # Output folder
# One or more GenomeRunner Web results data folders.
dname<-"data//more15_vs_tfbsEncode//"
mtx<-do.call("rbind", lapply(dname, function(fn) as.matrix(read.table(paste(fn, "matrix.txt", sep=""), sep="\t", header=T, row.names=1))))
# Optional: filter unused genomic features
# mtx<-mtx[grep("snp", rownames(mtx), ignore.case=T, invert=T), ]
mtx<-mtx.transform(mtx) # -log10 transform p-values
# Optional: adjust columns for multiple testing. See utils.R for the function definition.
# mtx<-mtx.adjust(mtx) 


## ----preprocessData, echo=TRUE, cache=TRUE, dependson='loadData'---------
dim(mtx) # Check original dimensions
# Define minimum number of times a row/col should have values above the cutoffs
numofsig<-1
cutoff<- -log10(0.01) # p-value significance cutoff
# What remains if we remove rows/cols with nothing significant
dim(mtx[apply(mtx, 1, function(x) sum(abs(x)>cutoff))>=numofsig, 
        apply(mtx, 2, function(x) sum(abs(x)>cutoff))>=numofsig])
# Trim the matrix
mtx<-mtx[apply(mtx, 1, function(x) sum(abs(x)>cutoff))>=numofsig, 
         apply(mtx, 2, function(x) sum(abs(x)>cutoff))>=numofsig]


## ----preprocessCorrel, echo=c(-4, -5), dependson='preprocessData'--------
# rcorr returns a list, [[1]] - correl coeffs, [[3]] - p-values. Type - pearson/spearman
mtx.cor<-rcorr(as.matrix(mtx), type="spearman")
# Optionally, try kendall correlation
# mtx.cor[[1]]<-cor(as.matrix(mtx), method="kendall")


## ----epigenomicVisualization, echo=c(-1:-6), fig.cap='Epigenomic similarity heatmap', fig.show='hold', fig.height=6.5----
par(oma=c(5,0,0,5), mar=c(10, 4.1, 4.1, 5)) # Adjust margins
color<-colorRampPalette(c("blue","yellow")) # Define color gradient
# Adjust clustering parameters.
# Distance: "euclidean", "maximum","manhattan" or "minkowski". Do not use "canberra" or "binary"
# Clustering: "ward", "single", "complete", "average", "mcquitty", "median" or "centroid"
dist.method<-"euclidean"  
hclust.method<-"ward"
h<-heatmap.2(as.matrix(mtx.cor[[1]]), trace="none", density.info="none", col=color,distfun=function(x){dist(x, method=dist.method)}, hclustfun=function(x){hclust(x, method=hclust.method)}, cexRow=0.5, cexCol=0.5)


## ----defineClusters, echo=c(-1:-3), results='markup', fig=TRUE-----------
par(oma=c(1, 0, 0, 0), mar=c(20, 4.1, 4.1,2.1), cex=0.5)
# Plot the dendrogram only, limit y axis. attr(h$colDendrogram, "height") has the maximum height of the dendrogram.
plot(h$colDendrogram, ylim=c(0, 15)) 
# Cut the dentrogram into separate clusters. Tweak the height
abline(h=4) # Visually evaluate the height where to cut
c<-cut(h$colDendrogram, h=4) 
# Check the number of clusters, and the number of members.
for (i in 1:length(c$lower)){
  cat(paste("Cluster", formatC(i, width=2, flag="0"), sep=""), "has", formatC(attr(c$lower[[i]], "members"), width=3), "members", "\n")
}
# Output the results into a file
unlink(paste(rname, "clustering.txt", sep=""))
for (i in 1:length(c$lower)){ 
  write.table(paste(i, t(labels(c$lower[[i]])), sep="\t"), paste(rname, "clustering.txt", sep=""), sep="\t", col.names=F, row.names=F, append=T)
}


## ----defineGroups, echo=TRUE, dependson='defineClusters'-----------------
eset.labels<-character() # Empty vector to hold cluster labels
eset.groups<-numeric() # Empty vector to hold cluster groups
# Set the minimum number of members to be considered for the differential analysis
minmembers<-9
for (i in 1:length(c$lower)) { # Go through each cluster
  # If the number of members is more than a minimum number of members
  if (attr(c$lower[[i]], "members") > minmembers) { 
    eset.labels<-append(eset.labels, labels(c$lower[[i]]))
    eset.groups<-append(eset.groups, rep(i, length(labels(c$lower[[i]]))))
  }
}


## ----limmaOnClusters, echo=TRUE, warning=FALSE, results='hide', dependson='defineGroups'----
eset<-new("ExpressionSet", exprs=as.matrix(mtx[, eset.labels]))
# Make model matrix
design<-model.matrix(~ 0+factor(eset.groups)) 
colnames(design)<-paste("c", unique(eset.groups), sep="")
# Create a square matrix of counts of DEGs
degs.matrix<-matrix(0, length(c$lower), length(c$lower))
colnames(degs.matrix)<-paste("c", seq(1,length(c$lower)), sep="")
rownames(degs.matrix)<-paste("c", seq(1, length(c$lower)), sep="") 
# Tweak p-value and log2 fold change cutoffs
cutoff.pval<-0.05
cutoff.lfc<-log2(2)
unlink(paste(rname, "degs.txt", sep=""))
for(i in colnames(design)){ 
  for(j in colnames(design)){
    # Test only unique pairs of clusters
    if (as.numeric(sub("c", "", i)) < as.numeric(sub("c", "", j))) {
      # Contrasts between two clusters
      contrast.matrix<-makeContrasts(contrasts=paste(i, j, sep="-"), levels=design)
      fit <- lmFit(eset, design) 
      fit2 <- contrasts.fit(fit, contrast.matrix)
      fit2 <- eBayes(fit2)
      degs<-topTable(fit2, number=dim(exprs(eset))[[1]], adjust.method="BH", p.value=cutoff.pval, lfc=cutoff.lfc)
      if(dim(degs)[[1]]>0) {
        print(paste(i, "vs.", j, ", number of degs:", dim(degs)[[1]]))
        # Keep the number of DEGs in the matrix
        degs.matrix[as.numeric(sub("c", "", i)), as.numeric(sub("c", "", j))]<-dim(degs)[[1]]
        # Average values in clusters i and j
        i.av<-rowMeans(matrix(exprs(eset)[rownames(degs), eset.groups == as.numeric(sub("c", "", i))], nrow=dim(degs)[[1]]))
        j.av<-rowMeans(matrix(exprs(eset)[rownames(degs), eset.groups == as.numeric(sub("c", "", j))], nrow=dim(degs)[[1]]))
        i.vs.j<-rep(paste(i,"vs.",j), dim(degs)[[1]])
        # Put it all together in a file, keeping columns with average transformed p-value being significant in at least one condition
        write.table(cbind(degs, i.vs.j, i.av, j.av)[abs(i.av) > -log10(cutoff.pval) || abs(j.av)> -log10(cutoff.pval),], paste(rname, "degs.txt", sep=""), sep="\t", col.names=NA, append=T)
      }
    }
  }
}


## ----maxminCorr, echo=c(-6), eval=TRUE, results='hide', dependson='preprocessCorrel'----
mtx.cor1<-mtx.cor[[1]]
# We don't need to consider perfect correlations, zero them out
diag(mtx.cor1)<-0
# Print top correlated parameters on screen
for (i in head(unique(mtx.cor1[order(mtx.cor1,decreasing=T)]))) {print(which(mtx.cor1 == i, arr.ind=T))}
unlink(paste(rname, "maxmin_correlations.txt", sep=""))
for (i in 1:ncol(mtx.cor1)) write.table(paste(colnames(mtx.cor1)[i],"correlates with",
                                              colnames(mtx.cor1)[which(mtx.cor1[i,] == max(mtx.cor1[i,]))], 
                                              "at corr. coeff.",formatC(mtx.cor1[i,which(mtx.cor1[i,] == max(mtx.cor1[i,]))]),
                                              "anticorrelates with",
                                              colnames(mtx.cor1)[which(mtx.cor1[i,] == min(mtx.cor1[i,]))],
                                              "at corr. coeff.",formatC(mtx.cor1[i,which(mtx.cor1[i,] == min(mtx.cor1[i,]))]),sep=","),
                                        paste(rname, "maxmin_correlations.csv", sep=""), append=T, sep=",", col.names=F, row.names=F) 


## ----enrichmentCutoffs, echo=c(-1), eval=TRUE, results='hide', fig.show='asis', dependson='preprocessCorrel'----
par(oma=c(1, 0, 0, 0), mar=c(5.1, 4.1, 4.1,2.1), cex=1)
# Define minimum number of times a row/col should have values above the cutoffs
numofsig<-1
dim(mtx) # Original dimensions
# Check summary and set p-value and variability cutoffs as 3rd quantiles of their distributions
summary(as.vector(abs(mtx))); cutoff.p<-summary(as.vector(abs(mtx)))[[5]]
summary(as.vector(apply(abs(mtx),1,sd))); cutoff.sd<-summary(as.vector(apply(abs(mtx),1,sd)))[[5]]
# Check visual distributions and set p-value and variability cutoffs manually
hist(as.vector(mtx), breaks=50, main="Distribution of -log10-transformed p-values", xlab="-log10-transformed p-values")
hist(c(as.vector(apply(mtx,1,sd)), as.vector(apply(mtx,2,sd))), breaks=50, main="Distribution of SD across rows and columns", xlab="SD")
cutoff.p<- -log10(0.05); cutoff.sd<-0.8


## ----trimCutoffs, echo=TRUE, results='hide', dependson='enrichmentCutoffs'----
mtx.gf<-mtx
# Remove rows/cols that do not show significant p-values and variability less than numofsig times
mtx.gf<-mtx.gf[apply(mtx.gf, 1, function(row){sum(abs(row)>cutoff.p)>=numofsig}), 
               apply(mtx.gf, 2, function(col){sum(abs(col)>cutoff.p)>=numofsig})]
mtx.gf<-mtx.gf[apply(mtx.gf, 1, sd)>cutoff.sd,
               apply(mtx.gf, 2, sd)>cutoff.sd]
dim(mtx.gf) # Dimensions after trimming


## ----enrichmentVisualization, echo=c(-1:-4, -7:-11, -13, -14), results='hide', warning=FALSE, fig=TRUE, fig.cap='Heatmap of the enrichment results', fig.show='asis', dependson='trimCutoffs', dependson='enrichmentCutoffs'----
# Adjust clustering parameters.
# Distance: "euclidean", "maximum","manhattan", "canberra", "binary" or "minkowski".
# Clustering: "ward", "single", "complete", "average", "mcquitty", "median" or "centroid"
dist.method<-"maximum"
hclust.method<-"ward"
# granularity=7
# my.breaks<-c(seq(min(mtx.gf), -cutoff.p, length.out=granularity), seq(cutoff.p, max(mtx.gf), length.out=granularity)) 
# my.breaks<-c(seq(min(mtx.gf), -2, length.out=granularity), seq(2, max(mtx.gf), length.out=granularity))
my.breaks<-c(seq(min(mtx.gf), max(mtx.gf)))
par(oma=c(10, 0, 0, 0), mar=c(5.1, 4.1, 4.1,2.1), cex=0.5)
color<-colorRampPalette(c("blue", "yellow"))
h<-heatmap.2(as.matrix(mtx.gf), distfun=function(x){dist(x,method=dist.method)}, hclustfun=function(x){hclust(x,method=hclust.method)}, dendrogram="none", breaks=my.breaks,  col=color, lwid=c(1.5,3), lhei=c(1.5,4), key=T,  symkey=T, keysize=0.01, density.info="none", trace="none",  cexCol=1.0, cexRow=0.8)


