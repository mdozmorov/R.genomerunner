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
    # from http://www.medcalc.org/calc/odds_ratio.php
    # Where zeros cause problems with computation of the odds ratio or its standard error, 0.5 is added to all cells (a, b, c, d) (Pagano & Gauvreau, 2000; Deeks & Higgins, 2010).
    tmp[(tmp$a == 0 | tmp$b == 0 | tmp$c == 0 | tmp$d == 0), ] <- tmp[(tmp$a == 0 | tmp$b == 0 | tmp$c == 0 | tmp$d == 0), ] + 0.5 # This makes calculations possible
    # create the new odds ratio column in the data frame
    mat$newOR = (tmp$a * tmp$d) / (tmp$b * tmp$c) 

    # from https://en.wikipedia.org/wiki/Fisher%27s_exact_test
    # asymptotic chi-square approximation is adequate if expected number of observations per cell is at least 5
    ind <- (tmp$a < 5 | tmp$b < 5 | tmp$c < 5 | tmp$d < 5) # Indexes where any of the cells is less than 5
    # Exclude unreliable calculations
    mat$newOR[ ind ] <- NA 
    
    # Now modify the pvalues based on this new odds ratio value.
    mat$newPval = mat$p_val
    mat$newPval[ mat$newPval < 1e-300 ] <- 1e-300 # Make 0 p-values to have some value, so they can be modified by OR
    mat$newPval[which(mat$newOR < 1)] = -1*mat$newPval[which(mat$newOR < 1)]
    # Exclude unreliable calculations
    mat$newPval[ ind ] <- NA 
    
    # Extract the new odds ratio and pvalues from the data frame and save to seperate files.
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
  write.table(as.data.frame(t(odds_mat)), file=paste(dirname(infile), "matrix_OR.txt", sep="/"), quote=FALSE, sep="\t")
  
  write.table(as.data.frame(t(pval_mat)), file=paste(dirname(infile), "matrix_PVAL.txt", sep="/"), quote=FALSE, sep="\t")
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
  system(paste("cat ",modfile," | sed -e 's/^-e //' > tmp && mv tmp ",modfile, sep=""))
  
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
    data[[as.character(parsed_GF$GF[i])]] <- v1[parsed_GF$dstart[i]:parsed_GF$dend[i],]
    
    tmp <- data.frame(a=as.numeric(paste(data[[i]][,"Observed"])), c=as.numeric(paste(data[[i]][,"Expected"])) )
    tmp$b = as.numeric(paste(parsed_GF$tot[i])) - tmp$a
    tmp$d = as.numeric(paste(parsed_GF$tot[i])) - tmp$c
    # from http://www.medcalc.org/calc/odds_ratio.php
    # Where zeros cause problems with computation of the odds ratio or its standard error, 0.5 is added to all cells (a, b, c, d) (Pagano & Gauvreau, 2000; Deeks & Higgins, 2010).
    tmp[(tmp$a == 0 | tmp$b == 0 | tmp$c == 0 | tmp$d == 0), ] <- tmp[(tmp$a == 0 | tmp$b == 0 | tmp$c == 0 | tmp$d == 0), ] + 0.5 # This makes calculations possible
    # create the new odds ratio column in the data frame
    data[[i]]$odds = (tmp$a * tmp$d) / (tmp$b * tmp$c) 
    
    # from https://en.wikipedia.org/wiki/Fisher%27s_exact_test
    # asymptotic chi-square approximation is adequate if expected number of observations per cell is at least 5
    ind <- (tmp$a < 5 | tmp$b < 5 | tmp$c < 5 | tmp$d < 5) # Indexes where any of the cells is less than 5
    # Exclude unreliable calculations
    data[[i]]$odds[ ind ] <- NA
  }
  
  # Now we have all info we need to extract columns into matricies!
  # get the odds ratio data
  oddsRatios <- as.data.frame(lapply(data, "[[", "odds"))
  row.names(oddsRatios) = pdf
  # write odds ratio data to file
  ## Now save odds matrix to file
  write.table(as.data.frame((oddsRatios)), file=paste(dirname(infile), "matrix_OR.txt", sep="/"), quote=FALSE, sep="\t")
  
  # Now get the pvalues.  They are modified to have a negative sign IF the odds ratio is less than 1.  
  # If odds ratio is NaN or Inf or -Inf or NA then nothing is done to the pvalues.  
  # this happens when the contingency table data contains zeros.
  modPvals <- lapply(data, "[[", "p.val")
  modPvals <- lapply(modPvals, "paste")
  modPvals <- as.matrix(as.data.frame(lapply(modPvals, "as.numeric")))
  modPvals[ modPvals < 1e-300 ] <- 1e-300 # Make 0 p-values to have some value, so they can be modified by OR
  row.names(modPvals) = pdf
  #modify pvalues to be negative if oddsRatios are less than 1.
  bool <- oddsRatios < 1
  bool[!is.finite(bool)] <- FALSE
  modPvals[bool] = -1*modPvals[bool]
  # Set unreliable p-values to NA
  modPvals[ is.na(oddsRatios) ] <- NA
  
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
