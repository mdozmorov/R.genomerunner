### Amy Olex
### 2/26/15
### Script to format GenomeRunner v2 files and extract out the odd ration column

#infile="GR_v2_log.txt"
#getV2OddsRatioMatrix(infile)
#currently this function assumes TCGA as it crops the row names.  I can take the cropping out if you want.
getV2OddsRatioMatrix <- function(infile){

  ## load file
  v2 <- read.delim(infile)
  
  ## Get GF names
  GF <- gsub("^###","",subset(v2,grepl("^###.*",foi_name))$foi_name)
  num_GF <- length(GF)
  ## Get the rest of the matrix without GF names
  mat <- subset(v2,!grepl("^###.*",foi_name))
  
  ## Get length of matrix
  mat_len <- dim(mat)[1]
  
  ## Get number of unique patient IDs
  num_foi <- length(unique(mat$foi_name))
  
  ## Ensure the number of rows in mat equals num_foi*num_GF
  ## This is needed to convert rows of data into matricies properly.
  if(num_foi*length(GF)==mat_len){
    ## Prior to converting to a matrix, re-calculate the odds ratios and correct for zeros
    tmp <- data.frame(a=as.numeric(paste(mat$foi_obs)))
    tmp$c = as.numeric(paste(mat$bg_obs)) - tmp$a
    tmp$b = as.numeric(paste(mat$n_fois)) - tmp$a
    tmp$d = as.numeric(paste(mat$n_bgs)) - as.numeric(paste(mat$n_fois)) - tmp$c
    ## Perform Fisher's exact test and store all the results
    tmp1 <- apply(tmp, 1, function(x) fisher.test(matrix(unlist(x), 2, 2)))
    # Keep regular OR
    mat$OR <- sapply(tmp1, function(x) {ifelse(x$conf.int[1] < 1 & x$conf.int[2] > 1, 1, x$estimate)})
    mat$OR[ is.infinite(mat$OR) ] <- .Machine$integer.max # Set infinite ORs to maximum number
    mat$OR[ mat$OR == 0 ] <- .Machine$double.eps # Set zero ORs to minimum number
    # If OR confidence interval includes 1, then OR is not significant, set to 1
    mat$newOR <- sapply(tmp1, function(x) {ifelse(x$conf.int[1] < 1 & x$conf.int[2] > 1, 1, ifelse(x$estimate < 1, x$conf.int[2], x$conf.int[1]))})
    mat$newOR[ is.infinite(mat$newOR) ] <- .Machine$integer.max # Set infinite ORs to maximum number
    mat$newOR[ mat$newOR == 0 ] <- .Machine$double.eps # Set zero ORs to minimum number
    
    # Now modify the pvalues based on this new odds ratio value.
    mat$newPval <- sapply(tmp1, "[[", "p.value") # Store p-values
    mat$newPval[ mat$newPval < 1e-300 ] <- 1e-300 # Make 0 p-values to have some value, so they can be modified by OR
    mat$newPval[which(mat$newOR < 1)] = -1*mat$newPval[which(mat$newOR < 1)]
    
    # Extract the new odds ratio and pvalues from the data frame and save to seperate files.
    oddsOR_mat=as.data.frame(matrix(mat$OR,nrow=num_foi,ncol=num_GF))
    names(oddsOR_mat) <- GF
    row.names(oddsOR_mat) <- mat$foi_name[1:num_foi] #substring(mat$foi_name[1:num_foi], 1, 15)
    
    odds_mat=as.data.frame(matrix(mat$newOR,nrow=num_foi,ncol=num_GF))
    names(odds_mat) <- GF
    row.names(odds_mat) <- mat$foi_name[1:num_foi] #substring(mat$foi_name[1:num_foi], 1, 15)
    
    pval_mat=as.data.frame(matrix(mat$newPval,nrow=num_foi,ncol=num_GF))
    names(pval_mat) <- GF
    row.names(pval_mat) <- mat$foi_name[1:num_foi]
  } else {
    print("Error: invalid number of foi or GF. num_foi*num_GF != num_mat_rows.")
    print(paste("num_foi=", num_foi, " num_GF=", length(GF), " num_mat_rows=", mat_len, sep=""))
    }
  
  ## Now save odds matrix and pval matrix to file
  ##write.table(t(odds_mat), file=paste(infile, ".OR", sep=""), quote=FALSE, sep="\t")
  write.table(as.data.frame(t(oddsOR_mat)), file=paste(dirname(infile), "matrix_ORFULL.txt", sep="/"), quote=FALSE, sep="\t")
  write.table(as.data.frame(t(odds_mat)), file=paste(dirname(infile), "matrix_OR.txt", sep="/"), quote=FALSE, sep="\t")
  write.table(as.data.frame(t(pval_mat)), file=paste(dirname(infile), "matrix_PVAL.txt", sep="/"), quote=FALSE, sep="\t")
  # Combine P-values and odds ratios
  comb_mat <- mtx.untransform(log2(odds_mat) * mtx.transform(pval_mat))
  write.table(as.data.frame(t(comb_mat)), file=paste(dirname(infile), "matrix_COMB.txt", sep="/"), quote=FALSE, sep="\t")
}


##########
## v1 file formatting
library(stringr)

## First!  Add the expected header line to the file.  I'm assuming these lines will NEVER change for the version 1 files.
## The data should have the header line 'foi_name    Observed  Expected  Diff	p-val	PCC	Obs/Tot'
## if it des not, this function will not work.

# infile="test_TFBS-DNAse-knownGene_LOG.gr.txt"
# infile="250bp_group1_f_TfbsClust3_LOG.gr.txt"
# infile="GR_v1_log"
# getV1OddsRatioPvalMatrix(infile)

getV1OddsRatioPvalMatrix <- function(infile){
  
  ## First we have to do some raw file formatting
  ## Create the command to insert a header line for the read.delim function
  modfile <- paste(infile, "_modified", sep="")
  cmd <- paste("echo -e \"foi_name\tObserved\tExpected\tDiff\tp-val\tPCC\tObs/Tot\" | cat - ",infile," > ",modfile, sep="\t")
  ## execute the command and then remove the extraneous "-e " in the file (not sure why it puts it there, doing it on terminal doesn't have this issue)
  system(cmd)
  system(paste("cat ",modfile," | sed -e 's/^-e //' > tmp.txt && mv tmp.txt ",modfile, sep=""))
  
  # Now, load file after manually copying the header line into the file:
  v1 <- read.delim(modfile)

  # Get a list of GFs:
  GF <- subset(v1,grepl("^Features analyzed:.*",foi_name))
  
  #parse out GF name and total number
  parsed_GF <- as.data.frame(str_match(GF$foi_name, "Features analyzed: (.*) \\(Total (.*)\\)"))[,2:3]
  names(parsed_GF) <- c("GF", "tot")
  
  #add in column for line number
  parsed_GF$line <- as.numeric(row.names(GF))
  
  #add data start line number
  parsed_GF$dstart <- parsed_GF$line+9
  
  if(dim(GF)[1] == 1){
    parsed_GF$dend <- dim(v1)[1]
  } else { 
    #get data end numbers (next line - 2)
    tmp <- parsed_GF$line[2:length(parsed_GF$line)]-2
    parsed_GF$dend <- c(tmp, dim(v1)[1])
  }
  
  #get patients/disease features
  pdf <- as.character(v1[parsed_GF$dstart[1]:parsed_GF$dend[1],]$foi_name)
  
  data <- list()
  
  #add data to the list and calculate odds ratio and insert column
  for(i in 1:dim(parsed_GF)[1]){
    data[[i]] <- v1[parsed_GF$dstart[i]:parsed_GF$dend[i],]
    
    tmp <- data.frame(a=as.numeric(paste(data[[i]][,"Observed"])), c=as.numeric(paste(data[[i]][,"Expected"])) )
    tmp$b = as.numeric(paste(parsed_GF$tot[i])) - tmp$a
    tmp$d = as.numeric(paste(parsed_GF$tot[i])) - tmp$c
    ## Perform Fisher's exact test and store all the results
    tmp1 <- apply(tmp, 1, function(x) fisher.test(matrix(unlist(x), 2, 2)))
    # If OR confidence interval includes 1, then OR is not significant, set to 1
    data[[i]]$newOR <- sapply(tmp1, function(x) {ifelse(x$conf.int[1] < 1 & x$conf.int[2] > 1, 1, ifelse(x$estimate < 1, x$conf.int[2], x$conf.int[1]))})
    data[[i]]$newOR[ is.infinite(data[[i]]$newOR) ] <- .Machine$integer.max # Set infinite ORs to maximum number
    data[[i]]$newOR[ data[[i]]$newOR == 0 ] <- .Machine$double.eps # Set zero ORs to minimum number
    
    # Now modify the pvalues based on this new odds ratio value.
    data[[i]]$newPval <- sapply(tmp1, "[[", "p.value") # Store p-values
    data[[i]]$newPval[ data[[i]]$newPval < 1e-300 ] <- 1e-300 # Make 0 p-values to have some value, so they can be modified by OR
    data[[i]]$newPval[which(data[[i]]$newOR < 1)] = -1*data[[i]]$newPval[which(data[[i]]$newOR < 1)]
  }
  
  # Now we have all info we need to extract columns into matricies!
  # get the odds ratio data
  oddsRatios <- as.data.frame(lapply(data, "[[", "newOR"))
  row.names(oddsRatios) = pdf
  colnames(oddsRatios) = parsed_GF$GF
  # write odds ratio data to file
  ## Now save odds matrix to file
  write.table(as.data.frame((oddsRatios)), file=paste(dirname(infile), "matrix_OR.txt", sep="/"), quote=FALSE, sep="\t")
  
  # Now get the pvalues.  They are modified to have a negative sign IF the odds ratio is less than 1.  
  # If odds ratio is NaN or Inf or -Inf or NA then nothing is done to the pvalues.  
  # this happens when the contingency table data contains zeros.
  modPvals <- lapply(data, "[[", "newPval")
  modPvals <- lapply(modPvals, "paste")
  modPvals <- as.matrix(as.data.frame(lapply(modPvals, "as.numeric")))
  row.names(modPvals) = pdf
  colnames(modPvals) = parsed_GF$GF
  # Now save the modified pvalues to a file
  write.table(as.data.frame((modPvals)), file=paste(dirname(infile), "matrix_PVAL.txt", sep="/"), quote=FALSE, sep="\t")
}

# ## Processing multiple folders
# 
# datadir <- "/Users/mikhail/Documents/Work/VCU_work/Edwin/meQTLs2/data.gr.filtered"
# 
# for (d in dir(datadir, pattern="^cr")) {
#   for (f in (list.files(paste(datadir, d, sep="/"), pattern="LOG.gr$", full.names=TRUE))) {
#     getV1OddsRatioPvalMatrix(f)
#     print(f)
#   }
# }
# 
# for (d in dir(datadir, pattern="^c")) {
#   for (f in (list.files(paste(datadir, d, sep="/"), pattern="^detailed", full.names=TRUE))) {
#     getV2OddsRatioMatrix(f)
#     print(f)
#   }
# }
