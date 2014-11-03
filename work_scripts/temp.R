# Not used
# Size vs. numofsigs analysis
length(list.files(pattern="^sel_")) # Hom many .gr files do we have?
tmp<-lapply(list.files(".", "^sel_"),function(i){
  m<-read.table(i,sep="\t",header=T,row.names=1) # Read in a current file as a matrix
  res <- t(as.data.frame(apply(m,2,function(x){sum(abs(x)>2)}))) # Get the number of significant elements columnwise
  rownames(res) <- i
  res
})
mg = function(x,y) {merge(x,y,all=T)} # Alias to merge function
tmp.merged<-Reduce(mg,tmp) # Apply it, as in http://stackoverflow.com/questions/8091303/merge-multiple-data-frames-in-a-list-simultaneously
tmp.merged # Row names will be messed up
tmp.m.summary<-t(as.matrix(apply(tmp.merged,2,function(i){sum(i,na.rm=T)}))) # Sum it up
tmp.m.size<-merge(gwassize,tmp.m.summary,all=T)
rownames(tmp.m.size)<-c("Size", "Num_Sig")
write.table(tmp.m.size,"clipboard",sep="\t")

# Mergingf

# Cory Reduce explanation
verbosePlus <- function(x,y) {
  print(paste(c(x,y)))
  x + y
}
Reduce(verbosePlus, 1:10)