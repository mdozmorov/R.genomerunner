## Convert a matrix of -log10-transformed p-values (with "-" indicating depletion) into raw linear scale
mtx.untransform<-function(x){
  if (nrow(x) == 0) return(x)
  tmp<- 1/10^abs(x) # -log10 transformation without sign
  for (i in 1:nrow(x)){
    for (j in 1:ncol(x)){
      if (x[i, j]<0) {tmp[i, j]<- -tmp[i,j]} # Add sign, if needed
    }
  }
  return(tmp)
}