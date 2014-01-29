# Convert transformed pvals into correlations
mtx<-mtx[apply(mtx, 1, function(row){sum(abs(row)>0)>=1}), ] #Remove rows with all zeros
mtx<-mtx[, apply(mtx, 2, function(col){sum(abs(col)>0)>=1})] #Remove columns with all zeros
mtx.cor<-rcorr(as.matrix(mtx), type="spearman") # [[1]] - correl coeffs, [[3]] - p-values. Type" pearson/spearman

# === Clustering of actual correlation coefficients
par(mar=c(10,6,6,5), oma=c(2,2,2,2)) # Make right and bottom margins larger
color<-colorRampPalette(c("blue","yellow"))
dist.method<-"euclidean"  # "euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski"
hclust.method<-"complete" # "ward", "single", "complete", "average", "mcquitty", "median" or "centroid"
#Maximum distance seems to best separate correlations. Centroid linkage further sharpen the separation
h<-heatmap.2(as.matrix(mtx.cor[[1]]), trace="none", density.info="none", col=color,distfun=function(x){dist(x, method=dist.method)}, hclustfun=function(x){hclust(x, method=hclust.method)}, cexRow=0.8, cexCol=0.8)

# Exploratory: Cllustering combinaations. Use to find bvisually best combinations of dist and hclust methods
dist.methods<-c("euclidean",  "manhattan", "minkowski", "canberra","maximum") # "binary",
hclust.methods<-c("ward", "single", "complete", "average", "mcquitty", "median", "centroid")
dev.off() # Clear graphic window
pdf(paste(dname, "cluster_combinations.pdf", sep=""))
# par(oma=c(5,0,0,5)) #Make right and bottom margins larger
for (d in dist.methods) {
  for (h in hclust.methods){
    # Correlations
    h<-heatmap.2(as.matrix(mtx.cor[[1]]), trace="none", density.info="none", col=color, distfun=function(x){dist(x, method=d)}, hclustfun=function(x){hclust(x, method=h)}, cexRow=0.5, cexCol=0.5, main=paste("Dist : ",d,"; Hclust : ",h))
  }
}
dev.off()

# Checking differenced among clusters
plot(h$rowDendrogram) # Just a dendrogram
# writeLines(colnames(mtx)[h$rowInd],"clipboard",sep="\t") # Order of  clustering
# writeLines(labels(rowDendrogram),"clipboard",sep="\t") # same thing
c<-cut(h$rowDendrogram,h=2.35) # Cut the dentrogram into separate cluster groups. Tweak the height
c
unlink(paste(dname, "clustering.txt", sep="")) # Clear file to write cluster groups
for (i in 1:length(c$lower)){ 
  write.table(paste(i, t(labels(c$lower[[i]])), sep="\t"), paste(dname, "clustering.txt", sep=""), sep="\t", col.names=F, row.names=F, append=T)
}

# Limma on cluster groups
eset.labels<-character() # Empty vector to hold cluster labels
eset.groups<-numeric() # Empty vector to hold cluster groups
for (i in 1:length(c$lower)) { # Go through each cluster
  if (attr(c$lower[[i]], "members") > 2) { # If the number of members is more than 2 (we need at least 3 members to test for statistically significant differences)
    eset.labels<-append(eset.labels, labels(c$lower[[i]])) # Append labels
    eset.groups<-append(eset.groups, rep(i, length(labels(c$lower[[i]])))) # Append numeric groups
  }
}
eset<-new("ExpressionSet", exprs=as.matrix(mtx[, eset.labels])) # Make eset out of eset.labels
design<-model.matrix(~ 0+factor(eset.groups)) # Make model matrix
colnames(design)<-paste("c", unique(eset.groups), sep="") # Colnames like "c1", "c2"
degs.matrix<-matrix(0, length(c$lower), length(c$lower)) # Square matrix of counts of DEGs. All clusters are present, even those not tested
colnames(degs.matrix)<-paste("c",seq(1,length(c$lower))); rownames(degs.matrix)<-paste("c",seq(1,length(c$lower))) # Name rows and cols by cluster name
cutoff.pval<-0.05 # p-value cutoff. Tweak
cutoff.lfc<-log(2) # log2 fold change cutoff. Tweak
for(i in colnames(design)){ # Loop through
  for(j in colnames(design)){ # each combination of clusters
    contrast.matrix<-makeContrasts(contrasts=paste(i, j, sep="-"), levels=design) # Contrasts between two clusters
    fit <- lmFit(eset, design) 
    fit2 <- contrasts.fit(fit, contrast.matrix)
    fit2 <- eBayes(fit2)
    degs<-topTable(fit2, number=dim(exprs(eset))[[1]], adjust.method="BH", p.value=cutoff.pval, lfc=cutoff.lfc)
    print(paste(i, "vs.", j, ", number of degs:", dim(degs)[[1]]))
    degs.matrix[as.numeric(sub("c","",i)), as.numeric(sub("c","",j))]<-dim(degs)[[1]] # Keep the number of DEGs in the matrix
    if(dim(degs)[[1]]>0) {
      i.av<-rowMeans(matrix(exprs(eset)[rownames(degs),eset.groups == as.numeric(sub("c","",i))], nrow=dim(degs)[[1]])) # Average value in cluster i
      j.av<-rowMeans(matrix(exprs(eset)[rownames(degs),eset.groups == as.numeric(sub("c","",j))], nrow=dim(degs)[[1]])) # Average value in cluster j
      i.vs.j<-rep(paste(i,"vs.",j), dim(degs)[[1]])
      write.table(cbind(degs, i.vs.j, i.av , j.av), paste(dname, "degs.txt", sep=""), sep="\t", col.names=NA, append=T)
    } # Write actual DEGs into a file
  }
}

# Checking max/min correlations
mtx.cor1<-mtx.cor[[1]]
diag(mtx.cor1)<-0 # We don't need to consider perfect correlations, zero them out
for (i in head(unique(mtx.cor1[order(mtx.cor1,decreasing=T)]))) print(which(mtx.cor1 == i, arr.ind=T)) # Print top correlated parameters
unlink(paste(dname, "maxmin_correlations.txt", sep="")) # Print best (anti)correlated values
for (i in 1:ncol(mtx.cor1)) write.table(paste(colnames(mtx.cor1)[i],"correlates with",
                                              colnames(mtx.cor1)[which(mtx.cor1[i,] == max(mtx.cor1[i,]))], 
                                              "at corr. coeff.",formatC(mtx.cor1[i,which(mtx.cor1[i,] == max(mtx.cor1[i,]))]),
                                              "anticorrelates with",
                                              colnames(mtx.cor1)[which(mtx.cor1[i,] == min(mtx.cor1[i,]))],
                                              "at corr. coeff.",formatC(mtx.cor1[i,which(mtx.cor1[i,] == min(mtx.cor1[i,]))]),sep="|"),
                                        paste(dname, "maxmin_correlations.txt", sep=""), append=T, sep="\t", col.names=F, row.names=F) 

# Clustering of the actual GFs
numofcols<-1 #Number of columns/rows with values above cutoff, if less, trim. Tweak it
dim(mtx) # Original dimensions
summary(as.vector(abs(mtx))); cutoff.p<-summary(as.vector(abs(mtx)))[[5]] # P-value cutoff as the 3rd quantile
summary(as.vector(apply(abs(mtx),1,sd))); cutoff.sd<-summary(as.vector(apply(abs(mtx),1,sd)))[[5]] # SD cutoff
hist(as.vector(mtx)); summary(as.vector(mtx)) # P-values distribution
hist(as.vector(apply(mtx,1,sd))); summary(as.vector(apply(mtx,1,sd))) # Variability across rows distribution
cutoff.p<-4; cutoff.sd<-0.5 # Cutoffs for p-value and SD. Tweak based on tests above
mtx.gf<-mtx # Temporary matrix for clustering
mtx.gf<-mtx.gf[apply(mtx.gf, 1, function(row){sum(abs(row)>cutoff.p)>=numofcols}), ] #Remove rows with less significant p-values
mtx.gf<-mtx.gf[apply(mtx.gf, 1, sd)>cutoff.sd, ] # Remove rows with less variability
mtx.gf<-mtx.gf[, apply(mtx.gf, 2, function(col){sum(abs(col)>cutoff)>=numofcols})] #Remove columts with no significant p-values
mtx.gf<-mtx.gf[, apply(mtx.gf, 2, sd)>cutoff.sd] # Remove cols with less variability
dim(mtx.gf) # Dims after trimming○
dist.method<-"maximum"  # "euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski"
hclust.method<-"ward" # "ward", "single", "complete", "average", "mcquitty", "median" or "centroid"
granularity=6
my.breaks<-c(seq(min(mtx.gf), -cutoff.p, length.out=granularity), seq(cutoff.p, max(mtx.gf), length.out=granularity)) #Breaks are going from min(mtx) to -cutoff, then black then from cutoff to max(mtx)
par(oma=c(10,0,0,5)) # Make right and bottom margins larger
color<-colorRampPalette(c("blue", "yellow"))
heatmap.2(as.matrix(mtx.gf), distfun=function(x){dist(x,method=dist.method)}, hclustfun=function(x){hclust(x,method=hclust.method)}, dendrogram="column", breaks=my.breaks,  col=color, lwid=c(1.5,3), lhei=c(1.5,4), key=T,  symkey=T, keysize=0.01, density.info="none", trace="none",  cexCol=1.0, cexRow=0.8)#, cellnote=formatC(mtx, digits=2), notecex=1, notecol='yellow') 
# No breaks for better contrast, but not to judge significance by color
heatmap.2(as.matrix(mtx.gf), distfun=function(x){dist(x,method=dist.method)}, hclustfun=function(x){hclust(x,method=hclust.method)}, dendrogram="none", col=color ,lwid=c(1.5,3), lhei=c(1.5,4), key=T,  symkey=T, keysize=0.01, density.info="none", trace="none",  cexCol=1.0, cexRow=0.8, font=2)#, cellnote=formatC(mtx, digits=2), notecex=1, notecol='yellow') 





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