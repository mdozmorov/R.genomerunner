# Convert a mantrix of raw p-values (with "-" indicating depletion) into -log10-transformed
mtx.transform<-function(x){
  tmp<- -log10(abs(x)) # -log10 transformation without sign
  for (i in 1:nrow(x)){
    for (j in 1:ncol(x)){
      if (x[i, j]<0) {tmp[i, j]<- -tmp[i,j]} # Add sign, if needed
    }
  }
  return(tmp)
}

# Correct for multiple testing the matrix of transformed p-values
mtx.adjust.log10<-function(x){ 
  tmp<- -log10(apply(1/10^abs(x), 2, p.adjust)) # Adjust absolute p-values for multiple testing and convert them back to -log10
  for (i in 1:nrow(x)){
    for (j in 1:ncol(x)){
      if (x[i,j]<0) {tmp[i,j]<- -tmp[i,j]} # Add sign, if neede
    }
  }
  return(tmp)
}

# Adjust raw p-values for mutliple testing
mtx.adjust.raw <- function(x, method="fdr") {
  tmp <- x; # Keep signs in the original version
  tmp.adj <- apply(abs(x), 2, p.adjust, method=method) 
  for (i in 1:nrow(x)){
    for (j in 1:ncol(x)){
      if (tmp[i, j] < 0) tmp.adj[i, j] <- -1*tmp.adj[i, j]
    }
  }
  return(tmp.adj)
}

barplot1<-function(mtx, bottom, l.x, l.y, ...){
  par(mar=c(bottom,5,2,2)+0.1)
  groupcolors<-c("blue","red") #c("yellow2","steelblue3","steelblue3","springgreen)
  b<-barplot(as.matrix(t(mtx)),beside=T,  ylab="-log10(p-value)\nnegative = underrepresentation",col=groupcolors,space=c(0.2,1),cex.names=1,las=2) #,names.arg=groupnames,legend.text=colnames(mtx),args.legend=list(x=7,y=4))
  lines(c(0,100),c(-log10(0.01),-log10(0.01)),type="l",lty="dashed",lwd=2)
  lines(c(0,100),c(log10(0.01),log10(0.01)),type="l",lty="dashed",lwd=2)
  
  legend(l.x,l.y,colnames(mtx),fill=groupcolors,cex=1)
  
}

