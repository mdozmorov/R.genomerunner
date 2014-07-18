# Online version. http://www.hiv.lanl.gov/content/sequence/HEATMAP/heatmap.html

library(gplots) # install.packages("gplot")
library(RColorBrewer) # of source("http:/bioconductor.org/biocLite.R") biocLite("RColorBrewer")
opar<-par(no.readonly=T) #Save original settings
par(oma=c(10,0,0,15)) #Make right and bottom margins larger
color<-rev(brewer.pal(11,'RdYlGn')) #Red-yellow-green gradient
color<-brewer.pal(6,'Reds') #Only for overrepresentation, for positive numbers, low=white, high=red intensity
color<-greenred #Standard green-black-red palette
# color<-colorRampPalette(c("red","white","green"))(10)

# === GenomeRunner correlations
mtx.raw <- as.data.frame(read.table("clipboard", sep="\t", header=T, row.names=1)) #Matrix of unadjusted p-values
mtx <- as.data.frame(read.table("clipboard", sep="\t", header=T, row.names=1)) #Matrix to work with, should be the same dimensions as mtx.raw
#omtx<-mtx #Keeps original working matrix, in case subsequent transformations would be unsatisfactory
mtx<-mtx.raw #Revert back changes to start over with

mtx.min<-min(mtx,na.rm=T) #Check Min/max, to know your limits. Do it after each transformation step.
min(mtx[mtx!=min(mtx)],na.rm=T) #Second to minimum
mtx.max<-max(mtx,na.rm=T) #Maximum
max(mtx[mtx!=max(mtx)],na.rm=T) #Second to maximum
print(c("Min: ", mtx.min, "  Max: ", mtx.max))
hist(as.vector(unlist(mtx)),breaks=50) # Hisogram of p-values distribution
# Filtering based on the number cutoff
cutoff<-2 # Cutoff is a number
mtx<-mtx[apply(mtx,1,function(row){sum(abs(row)>cutoff)>=1}),] #Remove rows with all values below cutoff
mtx<-mtx[,apply(mtx,2,function(col){sum(abs(col)>cutoff)>=1})] #Remove columns with all values below cutoff
# Optional: Set everything between -cutoff - cutoff to 0, based on unadjusted matrix.
# mtx[abs(mtx)<=cutoff]<-0
# Filtering based on the % cutoff
cutoff<-50 # Cutoff is a %
mtx<-mtx[apply(mtx,1,function(row){sum((row > mtx.max*cutoff/100 | row < mtx.min*cutoff/100))>=1}),] #Remove rows with all values below cutoff
mtx<-mtx[,apply(mtx,2,function(col){sum((col > mtx.max*cutoff/100 | col < mtx.min*cutoff/100))>=1})] #Remove columns with all values below cutoff
#Removing rows/cols with zeros may be sufficient, check mtx dimensions. Next is filtering by rows/cols containing at least numofcols values above cutoff
numofcols<-1 #Number of columns/rows with values above cutoff, if less, trim. Tweak it
mtx<-mtx[apply(mtx,1,function(row){sum(abs(row)>cutoff)>=numofcols}),] #Remove rows with no significant p-values
mtx<-mtx[,apply(mtx,2,function(col){sum(abs(col)>cutoff)>=numofcols})] #Remove colunts with no significant p-values
#If the results unsatisfactory, revert mtx<-omtx, and start over
write.table(as.data.frame(mtx),"clipboard",sep='\t') #Write trimmed matrix back to clipboard

dist.method<-"maximum"  # "euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski"
hclust.method<-"average" # "ward", "single", "complete", "average", "mcquitty", "median" or "centroid"
#With breaks, greenred colors.
granularity=6
my.breaks<-c(seq(-12,0,length.out=granularity),seq(0,12,length.out=granularity))
cutoff<-2 # Cutoff for calculating range of breaks
my.breaks<-c(seq(min(mtx),0,length.out=granularity),seq(0,max(mtx),length.out=granularity)) #Breaks are going from min(mtx) to 0, then black then from 0 to max(mtx)
my.breaks<-c(seq(min(mtx[mtx!=min(mtx)]),0,length.out=granularity),seq(0,max(mtx[mtx!=max(mtx)]),length.out=granularity)) #Breaks are going from second min(mtx) to -cutoff, then black then from cutoff to the second max(mtx)
my.breaks<-c(seq(min(mtx),-cutoff,length.out=granularity),seq(cutoff,max(mtx),length.out=granularity)) #Breaks are going from second min(mtx) to -cutoff, then black then from cutoff to the second max(mtx)
# my.breaks<-c(seq(min(mtx), mtx.min*cutoff/100, length.out=granularity), seq(mtx.max*cutoff/100, max(mtx), length.out=granularity)) #Breaks are going from min(mtx) to max(mtx) using % cutoff
heatmap.2(as.matrix(mtx), distfun=function(x){dist(x,method=dist.method)}, hclustfun=function(x){hclust(x,method=hclust.method)}, breaks=my.breaks, col=greenred(2*granularity-1),lwid=c(1.5,3), lhei=c(1.5,4), key=T,  symkey=T, keysize=0.1, density.info="none", trace="none",  cexCol=1, cexRow=1) #,labRow=groupnames) 
# Just gradient
heatmap.2(as.matrix(mtx), distfun=function(x){dist(x,method=dist.method)}, hclustfun=function(x){hclust(x,method=hclust.method)}, col=color,lwid=c(1.5,3), lhei=c(1.5,4), key=T,  symkey=T, keysize=0.1, density.info="none", trace="none",  cexCol=1.0, cexRow=1.0) 

# To have p-values overlayed on the cells: cellnote=formatC(1/10^mtx, format="e", digits=2)
heatmap.2(as.matrix(mtx), distfun=function(x){dist(x,method=dist.method)}, hclustfun=function(x){hclust(x,method=hclust.method)}, breaks=my.breaks, col=greenred(2*granularity-1),lwid=c(1.5,3), lhei=c(1.5,4), key=T,  symkey=T, keysize=0.1, density.info="none", trace="none",  cexCol=1.5, cexRow=1.5, cellnote=formatC(1/10^abs(mtx), format="e", digits=2), notecex=1.2, notecol='yellow')

# === Without breaks, brewer colors. Best for enrichment only
postscript("SLE_vs_TFs_Clustering.ps")
par(oma=c(5,0,0,5)) #Make right and bottom margins larger
heatmap.2(as.matrix((mtx)),trace="none",col=brewer.pal(6,'Reds'),distfun=function(x){dist(x,method=dist.method)}, hclustfun=function(x){hclust(x,method=hclust.method)}, density.info="none",cexCol=1,cexRow=1)#, cellnote=formatC(1/10^abs(mtx), format="e", digits=2), notecol='darkgreen')
dev.off()

# === No clustering, just heatmap
heatmap.2(as.matrix(mtx),Rowv=F,Colv=F,trace="none",col=color, density.info="none",cexCol=1,cexRow=1)#, cellnote=formatC(1/10^abs(mtx), format="e", digits=2), notecol='darkgreen')
# === Heatmap without clustering, for pre-defined rows and columns
heatmap.2(as.matrix(mtx),Rowv=FALSE,Colv=FALSE,col=brewer.pal(9,'Reds'),trace="none",density.info="none",cellnote=100*mtx,notecex=1.5,notecol='darkgreen')

# === Cycling through all combinations
mtx[mtx>0]<-1 #Set non-zero elements to 1, for presence/absence indication.
dist.methods<-c("euclidean",  "manhattan", "binary", "minkowski") #"canberra","maximum",
hclust.methods<-c("ward", "single", "complete", "average", "mcquitty", "median", "centroid")
dev.off() # Clear graphic window
pdf("test2.pdf")
par(oma=c(5,0,0,5)) #Make right and bottom margins larger
for (d in dist.methods) {
  for (h in hclust.methods){
    heatmap.2(as.matrix(mtx),trace="none",col=color,distfun=function(x){dist(x,method=d)}, hclustfun=function(x){hclust(x,method=h)}, density.info="none",cexCol=1,cexRow=0.8, notecex=1, main=paste("Dist : ",d,"; Hclust : ",h)) # cellnote=mtx.raw, 
  }
}
dev.off()



# === Getting clustering order http://stackoverflow.com/questions/5320814/order-of-rows-in-heatmap
mtx<-t(mtx)
mtx.dist<-dist(mtx,method=dist.method)
mtx.hclust<-hclust(mtx.dist,method=hclust.method)
mtx.dendrogram<-as.dendrogram(mtx.hclust)
Rowv<-rowMeans(mtx,na.rm=T)
dendrogram<-reorder(mtx.dendrogram,Rowv)
rowInd<-rev(order.dendrogram(dendrogram))
di<-dim(mtx)
nc<-di[2L]
nr<-di[1L]
colInd<-1L:nc
mtx.ordered<-mtx[rowInd,colInd]
write.table(mtx.ordered,"clipboard",sep='\t')

# === For cross-correlation expreiment only
mtx<-read.table("clipboard", header=T, sep='\t', row.names=1, as.is=T)
mtx.raw<-mtx
mtx<-mtx.raw
dim(na.omit(mtx))
mtx<-na.omit(mtx) # Remove NA rows
mtx<-mtx[apply(mtx,1,function(row){sum(abs(row)>0)>=1}),] #Remove rows with all zeros
mtx<-mtx[,apply(mtx,2,function(col){sum(abs(col)>0)>=1})] #Remove columns with all zeros
# mtx<-mtx[,!sapply(as.data.frame(mtx),sd)==0] #Remove columns with SD=0
# mtx<-mtx[!sapply(as.data.frame(t(mtx)),sd)==0,] # Remove rows with SD=0
mtx<-mtx[apply(mtx,1,sd)!=0,] # Remove rows with SD=0
mtx<-mtx[,apply(mtx,2,sd)!=0] # Remove cols with SD=0
# Convert measurements into correlations
mtx=cor(mtx,method="p") # "pearson" (default), "kendall", or "spearman"
diag(mtx)<-0
# for (i in 1:80) mtx[i,i]<-0 # zero out diagonal 1's
# for (i in head(unique(mtx[order(mtx,decreasing=T)]))) print(which(mtx == i, arr.ind=T)) # Print top correlated parameters
# Print best correlated values
# for (i in 1:14) print(paste("<",rownames(mtx)[i],">","correlates best with",
#                             "<",colnames(mtx)[which(mtx[i,] == max(mtx[i,]))],">",
#                             "at Pearson's",formatC(mtx[which(mtx[i,] == max(mtx[i,]))])))
# as.data.frame(t(apply(mtx,1,function(row) {x<-sort(row,decreasing=T)[1]; c(formatC(x),names(x))})))
# === Correlation mosaic. Method: "pearson", "kendall", "spearman"
dist.method<-"euclidean"  # "euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski"
hclust.method<-"ward" # "ward", "single", "complete", "average", "mcquitty", "median" or "centroid"
#Maximum distance seems to best separate correlations. Centroid linkage further sharpen the separation
color<-colorRampPalette(c("blue","yellow"))
color<-colorRampPalette(c("green","black","red"))
h<-heatmap.2(as.matrix(mtx),trace="none",col=color,distfun=function(x){dist(x,method=dist.method)}, hclustfun=function(x){hclust(x,method=hclust.method)})



plot(h$rowDendrogram) # Just a dendrogram
writeLines(rownames(mtx)[h$rowInd],"clipboard",sep="\t") # Order of  clustering
writeLines(labels(rowDendrogram),"clipboard",sep="\t") # same thing
c<-cut(h$rowDendrogram,h=300) # Cut the dentrogram
writeLines(labels(c$lower[[5]]), "clipboard",sep="\t") # Iterate to get names in subclusters
write.table(h$carpet,"F:/111.txt",sep="\t") # Clustered matrix


)mtx<-t(na.omit(t(mtx))) # Remove NA, if any are still there
cutoff<-0.5 #Everything below 10^cutoff is trimmed. 2 equal power of 10, so everything below 0.01 is not
numofcols<-2 #Number of columns/rows with values above cutoff, if less, trim. Tweak it
mtx<-mtx[apply(mtx,1,function(row){sum(abs(row)>cutoff)>=numofcols}),] #Remove rows with no significant p-values
mtx<-mtx[,apply(mtx,2,function(col){sum(abs(col)>cutoff)>=numofcols})] #Remove colunts with no significant p-values
dist.method<-"euclidean"  # "euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski"
hclust.method<-"ward" # "ward", "single", "complete", "average", "mcquitty", "median" or "centroid"
granularity=10
my.breaks<-c(seq(-1,0,length.out=granularity),seq(0,1,length.out=granularity)) #Breaks are going from min(mtx) to -cutoff, then black then from cutoff to max(mtx)
par(oma=c(5,0,0,5)) #Make right and bottom margins larger
# To have p-values overlayed on the cells: cellnote=formatC(mtx, digits=2)
heatmap.2(as.matrix(mtx), distfun=function(x){dist(x,method=dist.method)}, hclustfun=function(x){hclust(x,method=hclust.method)}, dendrogram="none", breaks=my.breaks, col=greenred(2*granularity-1), lwid=c(1.5,3), lhei=c(1.5,4), key=T,  symkey=T, keysize=0.01, density.info="none", trace="none",  cexCol=1.0, cexRow=0.8)#, cellnote=formatC(mtx, digits=2), notecex=1, notecol='yellow') 
# GrYlRd gradient
heatmap.2(as.matrix(mtx), distfun=function(x){dist(x,method=dist.method)}, hclustfun=function(x){hclust(x,method=hclust.method)}, dendrogram="none", breaks=my.breaks, col=rev(brewer.pal(2*granularity-1,'RdYlGn')),lwid=c(1.5,3), lhei=c(1.5,4), key=T,  symkey=T, keysize=0.01, density.info="none", trace="none",  cexCol=1.0, cexRow=1)#, cellnote=formatC(mtx, digits=2), notecex=1, notecol='yellow') 
# Standardizind the data
mtx<-scale(mtx,center=T,scale=T) # Standardize columns


# === Dendrogram
library(cluster)
IAC=mtx # cor(mtx,method="spearman") # method = c("pearson", "kendall", "spearman"))
# hist(IAC,sub=paste("Mean=",format(mean(IAC[upper.tri(IAC)]),digits=3)))
# mean(IAC)
# # Here we see that the mean IAC in the unnormalized dataset, with no outlier samples removed, is 0.905.  There is a long tail to the left of the distribution, indicating the presence of possible outliers.

# Performing hierachical clustering (average linkage) using 1-IAC as a distance metric
cluster1=hclust(as.dist(1-IAC),method="ward")
plot(cluster1,cex=0.7,labels=colnames(mtx),main="Naive B panel 2")

