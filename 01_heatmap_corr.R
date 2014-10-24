## ----preprocessCorrel, echo=c(-4, -5), dependson='preprocessData'--------
# rcorr returns a list, [[1]] - correl coeffs, [[3]] - p-values. Type - pearson/spearman
library(Hmisc)
mtx.cor<-rcorr(as.matrix(mtx), type="spearman")
# Optionally, try kendall correlation
# mtx.cor[[1]]<-cor(as.matrix(mtx), method="kendall")

## ----epigenomicVisualization, echo=c(-1:-6), fig.cap='Epigenomic similarity heatmap', fig.show='hold', fig.height=6.5----
library(gplots)
dev.off()
pdf(paste(rname, "mtx_clustered.pdf", sep=""))
par(oma=c(9,0,0,9), mar=c(12, 4.1, 4.1, 7)) # Adjust margins
color<-colorRampPalette(c("blue","yellow")) # Define color gradient
# Adjust clustering parameters.
# Distance: "euclidean", "maximum","manhattan" or "minkowski". Do not use "canberra" or "binary"
# Clustering: "ward", "single", "complete", "average", "mcquitty", "median" or "centroid"
dist.method<-"euclidean"  
hclust.method<-"ward.D2"
h<-heatmap.2(as.matrix(mtx.cor[[1]]), trace="none", density.info="none", col=color,distfun=function(x){dist(x, method=dist.method)}, hclustfun=function(x){hclust(x, method=hclust.method)}, cexRow=1, cexCol=1)
write.table(h$carpet, paste(rname, "mtx_clustered.txt", sep="/"), sep="\t", col.names=NA)
dev.off()

# Exploratory: Clustering combinaations. Use to find visually best combinations of dist and hclust methods
dist.methods<-c("euclidean",  "manhattan", "minkowski", "maximum") # "binary", "canberra",
hclust.methods<-"ward.D2" # c("ward", "single", "complete", "average", "mcquitty", "median", "centroid")
dev.off() # Clear graphic window
unlink(paste(rname, "cluster_combinations.pdf", sep=""))
pdf(paste(rname, "cluster_combinations.pdf", sep=""))
# par(oma=c(5,0,0,5)) #Make right and bottom margins larger
for (d in dist.methods) {
  for (h in hclust.methods){
    # Correlations
    h<-heatmap.2(as.matrix(mtx.cor[[1]]), trace="none", density.info="none", col=color, distfun=function(x){dist(x, method=d)}, hclustfun=function(x){hclust(x, method=h)}, cexRow=1, cexCol=1, main=paste("Dist : ",d,"; Hclust : ",h))
  }
}
dev.off()

## ----defineClusters, echo=c(-1:-3), results='markup', fig=TRUE-----------
par(oma=c(1, 0, 0, 0), mar=c(12, 4.1, 4.1,2.1), cex=1)
# Plot the dendrogram only, limit y axis. attr(h$colDendrogram, "height") has the maximum height of the dendrogram.
plot(h$colDendrogram) #, ylim=c(0, 15)) 
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




# Exploratory p-value distribution for each column
par(oma=c(0, 0, 0, 0) , mar=c(15, 6, 4, 2)+0.1) #Make right and bottom margins larger
boxplot(mtx[, order(apply(mtx, 2, max), decreasing=T)], cex.axis=0.8, las=2, ylab="-log10(p-value)\nNegative - underrepresentation", boxlwd=2, whisklwd=2, staplelwd=2)
title(xlab="p-value distribution", line=5)
abline(h=2,col="red",lwd=2) ; abline(h=-2, col="blue",lwd=2) # Significance cutoff

# Exploratory: Variability/SD (change) of the columns
mtx.var<-as.matrix(apply(mtx, 2, sd)) 
barplot(t(mtx.var[order(mtx.var, decreasing=T)]), names.arg=rownames(mtx.var)[order(mtx.var, decreasing=T)], cex.names=0.8, las=2, ylab="Standard deviation of the -log10(p-value)")
hist(t(mtx.var), n=50, cex=0.8) # Variability distribution
quantile(mtx.var)
cutoff<-quantile(mtx.var)[4] # 75% percentile
length(rownames(mtx.var)[mtx.var>cutoff]) # How many features have top 75% variability

# Exploratory: Check if variability correlate with sample size
gwassize<-read.table("data//gwassize.txt", sep="\t",header=T) # Features and the number of SNPs
setdiff(rownames(mtx.var), colnames(gwassize)) # Check how features in two matrixes overlap
size.var<-merge(t(gwassize), mtx.var, by="row.names") # Merge
rownames(size.var)<-size.var[,1] # Reassign row names
size.var<-size.var[,-1] # Remove first column with row names, after reassignment
colnames(size.var)<-c("snp","var") # Proper column names
# size.var<-size.var[size.var[,'var'] > 1,] # Optional: cut low variability snps, hope for better regression
attach(size.var)
plot(snp, var, pch=15, col="blue", ylim=c(0,13), xlab="gwas size", ylab="Standard deviation of the -log10(p-value)") # scatterplot one vs. the other
text(snp, var, labels=rownames(size.var), pos=3,cex=0.8) # Add text labels
res<-lm(var~snp) # Linear regression
abline(res, col="red", lwd=3) # Linear regression
lines(lowess(snp,var), col="blue") # lowes regression
rcorr(snp,var) # Pearsons corr, and p-value
detach(size.var)

# Exploratory: Check variability across GFs
gf.var<-as.matrix(apply(mtx,1,sd))
rownames(gf.var)<-rownames(mtx)
barplot(t(gf.var), cex.names=0.8, las=2)
hist(t(gf.var),n=50) # lots of GFs with var ~0.5

# Exploratory: P-value clustering
dist.method<-"euclidean"  # "euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski"
hclust.method<-"single" # "ward", "single", "complete", "average", "mcquitty", "median" or "centroid"
tmp<-as.matrix(-log10(mtx.cor[[3]])) # -log10 Matrix of p-values 
tmp.max<-max(tmp[!is.na(tmp) & !is.infinite(tmp)]) # Maximum non-NA and non-Inf value
diag(tmp)<-tmp.max + 2*.Machine[['double.xmin']] # Set NA values to just above max
tmp[is.infinite(tmp)]<-tmp.max + .Machine[['double.xmin']] # Set inf values to max, but less diagonal
par(mar=c(10,6,6,5),oma=c(2,2,2,2)) #Make right and bottom margins larger
color<-colorRampPalette(c("blue","yellow"))
h<-heatmap.2(tmp,trace="none",col=color,distfun=function(x){dist(x,method=dist.method)}, hclustfun=function(x){hclust(x,method=hclust.method)}, density.info="none",cexCol=0.8,cexRow=0.8) # ,cellnote=formatC(tmp, format="e", digits=2), notecol='darkgreen')

# Experimental: Clustering by sum of sizes. Outer will create one-dimensional array of sums, unlist will take it out of lizs, matrix will reformat as matrix
gwassize.2<-matrix(unlist(outer(1:length(gwassize),1:length(gwassize),FUN=function(i,j) gwassize[i]+gwassize[j])),nrow=length(gwassize),ncol=length(gwassize))
rownames(gwassize.2)<-colnames(gwassize.2)<-names(gwassize) # Reassign names
h<-heatmap.2(gwassize.2,trace="none",col=colorRampPalette(c("blue","white","red")),distfun=function(x){dist(x,method=dist.method)}, hclustfun=function(x){hclust(x,method=hclust.method)}, density.info="none",cexCol=0.8,cexRow=0.8 ,cellnote=formatC(tmp, format="e", digits=2), notecol='darkgreen')
write.table(mtx.cor[[3]][names(gwassize)[h$rowInd],names(gwassize)[h$colInd]],"clipboard-128",sep="\t") # Write p-values on the orded of clustering
write.table(gwassize[h$colInd],"clipboard",sep="\t") # And sizes
# Merge p-values and sizes
gwassize.3<-merge(as.matrix(t(gwassize)),mtx.cor[[3]],by="row.names") # Append gwassize
gwassize.3<-(merge(as.matrix(t(gwassize)),t(gwassize.3),by="row.names",all.y=T)) # transpose and append again
write.table(gwassize.3,"clipboard-128",sep="\t")

# Experimental: Comparing clusters
dev.off(); dev.new()
hc<-hclust(dist(as.matrix(mtx.cor[[1]]), method="euclidean"),method="average")
plot(hc, hang = -1)

cor(cophenetic(h$rowDendrogram),cophenetic(hc))
compareHclust <- function(hc1,hc2) {
  c1 <- as.matrix(cophenetic(hc1))
  c2 <- as.matrix(cophenetic(hc2))
  stopifnot(all(sort(rownames(c1))==sort(rownames(c2))))
  c2 <- c2[rownames(c1),rownames(c1)]
  list(r2=cor(as.dist(c1),as.dist(c2)),
       distance=diag(cor(c1,c2)))
}