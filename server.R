suppressMessages(source("utils2.R")) # See the required packages there
suppressMessages(source("episimilarity.R"))
(source("genomerunner_d3heatmap.R"))
library(d3heatmap)
library(dendextendRcpp) # required for extracting the height from the dendrogram
library(tools)
library(colorRamps)
library(shinyBS)
library(scales)
(source("functions/mtx.degfs.R"))

# # Lukas paths
results.dir <- "/home/lukas/db_2.00_06-10-2015/results/test2/"
gfAnnot <- read.table("/home/lukas/genome_runner/db/gf_descriptions.txt",sep="\t",header=T)
# # Mikhail paths
# gfAnnot <- read.table("/Users/mikhail/Documents/Work/GenomeRunner/genome_runner/db/gf_descriptions.txt", sep="\t",header=T)
# results.dir <- "/Users/mikhail/Documents/Work/GenomeRunner/R.GenomeRunner/data/test_30x5matrix_nonsig/"
# results.dir <- "/Users/mikhail/Documents/Work/GenomeRunner/Paper-Similarity/data_GWASdb2_manual/bed_selected/renamed/gappedPeak/"

genomerunner.mode <- F
coloring.num = 50
shinyServer(function(input, output,session) {
  
  
  # Parse the GET query string
  get.results.dir <- reactive({
    if (genomerunner.mode){
      query <- parseQueryString(session$clientData$url_search)
      return(paste(results.dir, query$id,"/",sep = ""))
    }else{
      return(results.dir)
    }
  })
  
  # enrichment ------------------------------------------------------------------
  get.matrix <- reactive({
    # populate the enrichment table combobox
    file.names.enrichment <- file_path_sans_ext(list.files(paste(get.results.dir(),"enrichment/",sep="")))
    mtx <- load_gr_data(paste(get.results.dir(), input$cmbMatrix,sep=""))
  })
  
  get.adjust.matrix <- reactive({
    mtx <- get.matrix()
    sign.mtx <- apply(mtx, 2, sign)
    adjust.rownames <- rownames(mtx)
    total.columns <- ncol(mtx)
    # Join with annotations
    mtx <- data.frame(GF=adjust.rownames, mtx) # Attach GF names
    mtx <- left_join(mtx, gfAnnot[, c("file_name", "cell")], by=c("GF" = "file_name"))
    row.names(mtx) <-  adjust.rownames
    class(mtx$cell) <- "character"
    mtx$cell[ is.na(mtx$cell) ] <- "dummy_cell" # If some file names is not in the gfAnnot dataframe (e.g., user-provided data), 'cell' column will contain NAs. replace them with dummy text to allow FDR correction
    unique.cells <- unique(mtx$cell) # Keep unique cell types
    
    # Adjust for multiple testing and -log10 transform, if needed
    for (i in 1:total.columns) { # Process each column
      for (u.c in unique.cells) { # Adjust for multiple testing on per-cell-type basis
        # If the cell-specific subset have >1 row, perform correction for multiple testing
        if(sum(mtx$cell == u.c) > 1) {
          # i+1 because we added the GF column
          mtx[mtx$cell == u.c,i + 1] <- mtx.transform(apply(mtx.untransform(mtx[mtx$cell == u.c,i+1,drop=F]), 2, function(x) p.adjust(abs(x), method = input$cmbPvalAdjustMethod)))
        }
      }  
    }
  
    mtx <- sign.mtx*mtx[,2:total.columns+1]
  })
  
  output$heatmapEnrich <- renderD3heatmap({
    if (input$cmbMatrix == "matrix_PVAL.txt"){
      mtx <- get.adjust.matrix()
    }else{
      mtx <- get.matrix()
    }
    n_limit = 20
    # if n > 100, calculate SD for each row.
    if (nrow(mtx) > n_limit){
      #  calculate SD for each row.
      mtx.sd <- apply(mtx,1,sd)
      mtx.sd <- data.frame(mtx.sd)
      # sort rows based on SD
      mtx.sd.order <- mtx[order(mtx.sd,decreasing = T),]
      mtx.sd.order <- mtx.sd.order[1:n_limit,]
      mtx <- mtx.sd.order
    }
    
    
    coloring<-colorRampPalette(c("blue", "yellow", "red"))
    d3heatmap::d3heatmap(as.matrix(mtx),heatmap_options = list(hclust=function(tmp) {hclust(tmp, method = input$cmbClustMethod)}), colors = coloring(coloring.num), show_tip=FALSE,dendro.rds.path=paste(get.results.dir(),"heatmap.dend.rds", sep=""))
    
  })
  
  
  
  output$legendEnrich <- renderPlot({
    if (input$cmbMatrix == "matrix_PVAL.txt"){
      mtx <- get.adjust.matrix()
    }else{
      mtx <- get.matrix()
    }
    coloring<-colorRampPalette(c("blue", "yellow", "red"))
    my.breaks <- c(seq(min(mtx), 0, length.out=10), 0, seq(0, max(mtx), length.out=10)) # Breaks symmetric around 0
    
    plot(my.breaks,rep(1,length(my.breaks)),col=coloring(length(my.breaks)),pch=15,cex=10,main="Depletion/Enrichment Significance",ylab="",xlab="",yaxt="n",xaxt='n')
    if (input$cmbMatrix == "matrix_PVAL.txt"){
      axis(side = 1, at =  my.breaks, labels=scientific_format(2)(1/10^abs( my.breaks)),las=1)
    }else{
      axis(side = 1, at =  my.breaks, labels=scientific_format(2)(2^my.breaks),las=1)
    }
  })
  
  # generate the enrichment table and appends the gf.name information columns
  get.enrichment.table <- reactive({
    mtx <- read.csv(paste(get.results.dir(),input$cmbMatrix,sep = ""),sep="\t")
    selectedFOI <- 1
    selectedFOI <-input$cmbFOI
    if (input$cmbMatrix == "matrix_PVAL.txt"){
      mtx.adjust <- apply(mtx[selectedFOI], 2, function(x) p.adjust(abs(x), method = input$cmbPvalAdjustMethod))
      pval.sig <- rep("Not significant",nrow(mtx.adjust))
      pval.sig[ mtx[, selectedFOI] < 0.05 & mtx[, selectedFOI] > 0] <- "Overrepresented" 
      pval.sig[mtx[, selectedFOI] > -0.05 & mtx[, selectedFOI] < 0] <- "Underrepresented"
      
      mtx.table <- cbind(abs(mtx[selectedFOI]), 
                         pval.sig, 
                         mtx.adjust)
      colnames(mtx.table) <- c("P.value","Direction","P.adj")
    }else{ 
      # for odds ratio
      or.sig <- rep("Not significant", nrow(mtx))
      or.sig[mtx[, selectedFOI] > 1] <- "Enriched"
      or.sig[mtx[, selectedFOI] < 1] <- "Depleted"
      mtx.table <- cbind(abs(mtx[selectedFOI]),
                         or.sig)
      colnames(mtx.table) <- c("Odds.ratio","Direction")
    }
    mtx.table = cbind(GF.Name = rownames(mtx.table),mtx.table) # make the GF.name a column instead of just the rowname so we can left_join
    rownames(mtx.table) <- NULL
    # gfAnnot is loaded in utils2
    
    mtx.table <- left_join(mtx.table,gfAnnot[,(names(gfAnnot) %in% c("file_name", "cell", "cell_desc", "factor", "factor_desc", "source", "source_desc"))],
                           by=c("GF.Name"="file_name"))
    
    return(mtx.table)
  })
  
  output$tblEnrichment <-renderDataTable({
    num.char <- 50
    table.enrich <- get.enrichment.table()
    table.enrich <- apply(table.enrich,c(1,2),function(x) {
      if (!is.na(x) & nchar(x)>num.char){
        return(paste(substring(x,1,num.char),  "..."))
      } else{
        return(x)
      }
    })
  }, options = list( lengthMenu = list(c(10, 50, 100,-1), c('10', '50','100', 'All')),
                     pageLength = 50))
  
  
  
  output$downloadEnrichTable <- downloadHandler(
    filename = function() { 
      return("Enrichment_table.txt")
    },
    content = function(file) {
      write.table(x = get.enrichment.table(),file =  file ,sep = "\t",quote = F,row.names = F)
    }
  )
  
  get.annotation.table <- reactive({
    validate(need(try(mtx <- read.table(paste(get.results.dir(),"annotations/", input$cmbAnnotation,sep=""),header=T)),"Annotation results not available."))
    mtx
  })
  
  output$tblAnnotation <- renderDataTable({
    get.annotation.table()
  })
  
  outputOptions(output, "downloadEnrichTable", suspendWhenHidden=FALSE)
  
  ## enrichment up and down plots for single column --------------------------------------------------------------
  get.barplot.matrix <- reactive({
    # populate the enrichment table combobox
    file.names.enrichment <- file_path_sans_ext(list.files(paste(get.results.dir(),"enrichment/",sep="")))
    mtx <- load_gr_data(paste(get.results.dir(), input$cmbMatrix,sep=""))
    
  })
  
  
  getEnrichmentUp <- reactive({
    mtx <-  data.frame( get.barplot.matrix())
   
    
    selectedFOI <- 1
    selectedFOI <-input$cmbFOI
    updown.split =  0
    # check if any results to subset
    if (nrow(mtx)==0){
      return(mtx)
    }
    mtx.up <- subset(mtx[selectedFOI],mtx[selectedFOI] > updown.split)
  
    if (nrow(mtx.up)==0){
      return(mtx.up)
    }
    
    if (input$cmbMatrix == "matrix_PVAL.txt"){
      
      # do the pvalue adjustment 
      mtx.up.adjust <- apply(1/10^(abs(mtx.up)), 2, function(x) {p.adjust(x,  method = input$cmbPvalAdjustMethod)})
      mtx.up.adjust <- as.matrix(mtx.up.adjust)
      rownames(mtx.up.adjust) <- rownames(mtx.up); colnames(mtx.up.adjust) <- colnames(mtx.up);
      mtx.up <- mtx.transform(mtx.up.adjust);
    }
    
    
    # sort the results
    mtx.up.sorted <- mtx.up[order(mtx.up,decreasing = T),,drop=FALSE]
  })
  
  getEnrichmentDown <- reactive({
    mtx <-  data.frame( get.barplot.matrix())
    selectedFOI <- 1
    selectedFOI <-input$cmbFOI
    updown.split = 0
    
    if (nrow(mtx)==0){
      return(mtx)
    }
    mtx.down <- subset(mtx[selectedFOI],mtx[selectedFOI] < updown.split,drop = F)
    
    if (nrow(mtx.down)==0){
      return(mtx.down)
    }
    
    if (input$cmbMatrix == "matrix_PVAL.txt"){
      
      # do the pvalue adjustment 
      mtx.down.adjust <- apply(1/10^(abs(mtx.down)), 2, function(x) {p.adjust(x,  method = input$cmbPvalAdjustMethod)})
      mtx.down.adjust <- -as.matrix(mtx.down.adjust) # apply negative bc everything is downregulated here
      rownames(mtx.down.adjust) <- rownames(mtx.down); colnames(mtx.down.adjust) <- colnames(mtx.down);
      mtx.down <- mtx.transform(mtx.down.adjust);
    }
    mtx.down.sorted <- mtx.down[order(mtx.down,decreasing = F),,drop=FALSE]
  })
  
  # the same barplot is used for single FOI results and multiple FOI results
  output$pltEnrichUp <- renderPlot({
    mtx.up.sorted = getEnrichmentUp()
    if (nrow(mtx.up.sorted)==0){
      # plot raw text
      par(mar = c(0,0,0,0))
      plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
      text(x = 0.5, y = 0.5, paste("Nothing overrepresented to plot."),  cex = 1.6, col = "black")
      box()
      return()
    }
    par(mar = c(10,7,4.1,2.1))
    if (nrow(mtx.up.sorted) < 3){
     x.loc = barplot(as.matrix(t(head(mtx.up.sorted,30))), beside=T,col = "red3",
              space=c(0.2,1), cex.names=1, las=2, xaxt='n',main="Enriched epigenomic associations",xlim=c(0,10*nrow(mtx.up.sorted)),axes=F)
    } else{
     x.loc =  barplot(as.matrix(t(head(mtx.up.sorted,30))), beside=T,col = "red3",
              space=c(0.2,1), cex.names=1, las=2, xaxt='n',main="Enriched epigenomic associations",axes=F)
    }
    abline(a=0,b=0)   
    axis.values = seq(0,max(mtx.up.sorted),length.out = 10)
    # draw y labels
    if(input$cmbMatrix == "matrix_PVAL.txt"){
      mtext('P-value',side=2,line=4)
      axis(side = 2, at = axis.values,labels=scientific_format(2)(1/10^abs(axis.values)),las=1)
    } else{
      mtext('Odds-ratio',side=2,line=5)
      axis(side = 2, at = axis.values, labels=scientific_format(2)(2^axis.values),las=1)
    }
    # draw rotated x-labels
    axis(side=1, labels = FALSE,tick = F)
    text(x.loc, par("usr")[3], adj=c(1,1),srt = 45,
         labels = head(rownames(mtx.up.sorted),30), xpd = TRUE) 
  })
  
  output$pltEnrichDown <- renderPlot({
    mtx.down.sorted <- abs(getEnrichmentDown())
    if (nrow(mtx.down.sorted) == 0){
      # plot raw text
      par(mar = c(0,0,0,0))
      plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
      text(x = 0.5, y = 0.5, paste("Nothing underrepresented to plot."),  cex = 1.6, col = "black")
      box()
      return(mtx.down.sorted)
    }
    par(mar = c(10,7,4.1,2.1))
    if (nrow(mtx.down.sorted) < 3){
      x.loc = barplot(as.matrix(t(head(mtx.down.sorted,30))), beside=T,col = "green4",
              space=c(0.2,1), cex.names=1, las=2, names.arg=head(rownames(mtx.down.sorted),30),xaxt='n',main = "Depleted epigenomic associations",xlim=c(0,10*nrow(mtx.down.sorted)),axes = F)
    }else{
      x.loc = barplot(as.matrix(t(head(mtx.down.sorted,30))), beside=T,col = "green4",
              space=c(0.2,1), cex.names=1, las=2, names.arg=head(rownames(mtx.down.sorted),30),xaxt='n',main = "Depleted epigenomic associations", axes = F)
    }
    abline(a=0,b=0) 
    # draw y-labels
    axis.values = seq(0,max(mtx.down.sorted),length.out = 10)
    if(input$cmbMatrix == "matrix_PVAL.txt"){
      mtext('P-value',side=2,line=4)
      axis(side = 2, at = axis.values,labels=scientific_format(2)(1/10^abs(axis.values)), las=1)
    }
    else{
      mtext('Odds-ratio',side=2,line=5)
      axis(side = 2, at = axis.values, label=scientific_format(2)(2^(-axis.values)), las=1)
    }
    # draw rotated x-labels
    axis(side=1, labels = FALSE,tick = F)
    text(x.loc, par("usr")[3], srt = 45, adj=c(1,1),
         labels = head(rownames(mtx.down.sorted),30), xpd = TRUE) 
  })
  
  # episimilarity ---------------------------------------------------------------
  get.corr.matrix <- reactive({
    mtx <- load_gr_data(paste(get.results.dir(), input$cmbMatrix,sep=""))
    # If there are columns with SD=0, add jitter to it. Needed for pair-wise column correlation analysis (epigenomic similarity analysis). Only valid if there's more than 1 row
    if (nrow(mtx) > 1) {
      ind <- apply(mtx, 2, function(x) sd(x, na.rm=TRUE)) == 0 # Indexes of such columns
      if (sum(ind) > 0) {
        set.seed(1)
        mtx[, ind] <- jitter(mtx[, ind, drop=FALSE], factor=0.1)
      }
    }
    mtx <- scale(mtx)
    
   
    rcorr(as.matrix(mtx), type=input$cmbEpisimCorType)[[1]]
  })
  
  get.cor.hclust.dendrogram <- reactive({
    cor.mat <- get.corr.matrix() 
    as.dendrogram(hclust(as.dist(1-cor.mat), method=input$cmbClustMethod))
  })
  
  output$heatmapEpisim <- renderD3heatmap({ 
    mat <- get.matrix()
    validate(need(ncol(mat)>2,"Need at least 3 SNPs of interest files to perform clustering."))
    validate(need(nrow(mat)>4,"Need at least 5 genome features to perform clustering."))
    cor.mat <- get.corr.matrix()
    hclustergram <- get.cor.hclust.dendrogram()
    coloring<-colorRampPalette(c("blue", "yellow", "red"))
    d3heatmap::d3heatmap(as.matrix(cor.mat),heatmap_options = list(Rowv=hclustergram,Colv=hclustergram,keep.dendro=TRUE),colors = coloring(coloring.num),dendro.rds.path=paste(get.results.dir(),"heatmap.dend.rds", sep=""))
  })
  
  output$legendEpisim <- renderPlot({
    mat <- get.matrix()
    validate(need(ncol(mat)>2,""))
    validate(need(nrow(mat)>4,""))
    mtx <- get.corr.matrix()
    coloring<-colorRampPalette(c("blue", "yellow", "red"))
    color.range = seq(min(mtx),max(mtx), (max(mtx)-min(mtx))/coloring.num )
    plot(color.range,rep(1,coloring.num+1),col=coloring(coloring.num+1),pch=15,cex=10,main="Correlation Coefficient",ylab="",xlab="",yaxt="n")
  })
  
  # this function is cut out from the tblEpigenetics renderer. It is a long calculation that is only run when # of clusters changes
  calculate.clust <- reactive({
    cor.mat <- get.corr.matrix() # this line ensure that dendrogram is redrawn when heatmap is
    hclustergram <- get.cor.hclust.dendrogram() # ensures that dendrogram is redrawn when hclust method is changed
    
    dend = readRDS(file = paste(get.results.dir(), "heatmap.dend.rds",sep=""))
    cl_num <- input$sldEpisimNumClust # Empirically set desired numter of clusters
    hcut <- heights_per_k.dendrogram(dend)[cl_num] # extract the height to cut based on # of groups
    # get the cluster labels
    mtx.clust <- dend %>% mtx.clusters(height=hcut, minmembers=3)
    mtx = load_gr_data(paste(get.results.dir(), input$cmbMatrix,sep="")) # load the original matrix
    is.OR = T
    if(input$cmbMatrix == "matrix_PVAL.txt"){
      is.OR = F
    }
    mtx.deg <- suppressWarnings(mtx.degfs(mtx[, mtx.clust$eset.labels], mtx.clust, label="broadPeak2",isOR = is.OR))
    
    updateSelectInput(session,"cmbEpigenetics","Select which epigenetic table to render",choices = names(mtx.deg))
    mtx.deg.path = paste(get.results.dir(),"mtx.deg.episim.RDS",sep = "")
    saveRDS(mtx.deg,file=mtx.deg.path)
    return(mtx.deg)
  })
  
  get.epigenetics.table <- reactive({
    mtx.deg <- calculate.clust()
    # check if any results were returned
    if (is.null(names(mtx.deg))){ 
      return(data.frame(NoResult="There is nothing signficant to show"))
    }
    selectedCor = names(mtx.deg)[1]
    # save last selected value
    if (input$cmbEpigenetics != "Results not ready yet.") {
      if (input$cmbEpigenetics %in% names(mtx.deg)){
        selectedCor = input$cmbEpigenetics
      }
      updateSelectInput(session,"cmbEpigenetics",choices=names(mtx.deg),selected = selectedCor)
    }
    #convert values to numeric form for sorting purposes
    for(x in list("adj.p.val",3,4)){
      mtx.deg[[selectedCor]][[x]] <- as.numeric(mtx.deg[[selectedCor]][[x]])
    }
    mtx.deg[[selectedCor]]
  })
  
  output$tblEpigenetics <-renderDataTable({
    mtx <- get.matrix()
    validate(need(ncol(mtx)>2,"Need at least 3 SNPs of interest files to perform clustering."))
    validate(need(nrow(mtx)>4,"Need at 5 least genome features to perform clustering."))
    validate(need(try(get.epigenetics.table()),"Try a different clustering method."))
    get.epigenetics.table()
  },options = list( lengthMenu = list(c(10, 50, 100,-1), c('10', '50','100', 'All')),
                    pageLength = 50))
  
  output$downloadEpigenetics <- downloadHandler(
    filename = function() { 
      return("Epigenetics_table.txt")
    },
    content = function(file) {
      write.table(x = get.epigenetics.table(),file =  file ,sep = "\t",quote = F,row.names = F)
    }
  )
  
  output$pltDend <- renderPlot({ 
    mtx <- get.matrix()
    validate(need(ncol(mtx)>2,""))
    validate(need(nrow(mtx)>4,""))
    
    cor.mat <- get.corr.matrix() # this line ensure that dendrogram is redrawn when heatmap is
    hclustergram <- get.cor.hclust.dendrogram() # ensures that dendrogram is redrawn when hclust method is changed
    
    dend = readRDS(file = paste(get.results.dir(), "heatmap.dend.rds",sep=""))
    plot(as.dendrogram(dend, hang=-1)) # Plot dendrogram
    cl_num <- input$sldEpisimNumClust # Empirically set desired numter of clusters
    cols <- rainbow(cl_num) # Make colos for clusters
    print(cl_num)
    hcut <- heights_per_k.dendrogram(dend)[cl_num] # extract the height to cut based on # of groups
    # get the cluster labels
    mtx.clust <-validate(need(try(dend %>% mtx.clusters(height=hcut, minmembers=3)),"Try using a lower number of clusters"))
    # write.table(as.data.frame(mtx.clust), "/home/lukas/clustering_all.txt", sep="\t", row.names=FALSE, quote=FALSE)
    mtx = load_gr_data(paste(get.results.dir(), input$cmbMatrix,sep="")) # load the original matrix
    
    # Define the clusters by rectangles
    validate(need(try(rect.hclust( as.hclust(dend), k=cl_num, border=cols)),"Try a different clustering method."))
    
  })
  
  
  # --download button code --------------------------------------------------
  
  
  output$downloadEnrichHeatmap <- downloadHandler(
    filename = function() { 
      return("EnrichmentHeatmap.pdf")
    },
    content = function(file) {
      pdf(file=file,width=10,height=10)
      par(mfrow=c(3,3))
      if (input$cmbMatrix == "matrix_PVAL.txt"){
        mtx <- get.adjust.matrix()
      }else{
        mtx <- get.matrix()
      }
      n_limit = 20
      # if n > 100, calculate SD for each row.
      if (nrow(mtx) > n_limit){
        #  calculate SD for each row.
        mtx.sd <- apply(mtx,1,sd)
        mtx.sd <- data.frame(mtx.sd)
        # sort rows based on SD
        mtx.sd.order <- mtx[order(mtx.sd,decreasing = T),]
        mtx.sd.order <- mtx.sd.order[1:n_limit,]
        mtx <- mtx.sd.order
      }
     
      heatmap.2(as.matrix(mtx),hclust=function(tmp) {hclust(tmp, method = input$cmbClustMethod)},density.info = 'none',main = "Enrichment Heatmap",
                margins = c(10,10),srtRow = -45,srtCol = 45)
      dev.off()
    },
    contentType = 'application/pdf'
  )
  
  output$downloadEpisimHeatmap <- downloadHandler(
    filename = function() { 
      return("EpisimilarityHeatmap.pdf")
    },
    content = function(file) {
      pdf(file=file,width=10,height=10)
      mat <- get.matrix()
      validate(need(ncol(mat)>2,"Need at least 3 SNPs of interest files to perform clustering."))
      validate(need(nrow(mat)>4,"Need at least 5 genome features to perform clustering."))
      cor.mat <- get.corr.matrix()
      hclustergram <- get.cor.hclust.dendrogram()
      coloring<-colorRampPalette(c("blue", "yellow", "red"))
      heatmap.2(as.matrix(cor.mat),Colv = hclustergram,margins = c(10,10),density.info = 'none',srtRow = -45,srtCol = 45,main = "Episimilarity Heatmap")
      dev.off()
    },
    contentType = 'application/pdf'
  )
  
  output$downloadEnrichBarPDF <- downloadHandler(
    filename = function() { 
      return("EnrichmentBarPlot.pdf")
    },
    content = function(file) {
      pdf(file=file,width=9,height=5)
      # print EnrichmentUp to PDF
      mtx.up.sorted = getEnrichmentUp()
      if (nrow(mtx.up.sorted)==0){
        # plot raw text
        par(mar = c(0,0,0,0))
        plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
        text(x = 0.5, y = 0.5, paste("Nothing overrepresented to plot."),  cex = 1.6, col = "black")
        box()
      }else{
        par(mar = c(10,7,4.1,2.1))
        if (nrow(mtx.up.sorted) < 3){
          x.loc <-barplot(as.matrix(t(head(mtx.up.sorted,30))), beside=T,col = "red3",
                  space=c(0.2,1), cex.names=1, las=2, names.arg=head(rownames(mtx.up.sorted),30),xaxt='n',main="Enriched epigenomic associations",xlim=c(0,10*nrow(mtx.up.sorted)),axes=F)
        } else{
          x.loc <-barplot(as.matrix(t(head(mtx.up.sorted,30))), beside=T,col = "red3",
                  space=c(0.2,1), cex.names=1, las=2, names.arg=head(rownames(mtx.up.sorted),30),xaxt='n',main="Enriched epigenomic associations",axes=F)
        }
        abline(a=0,b=0)   
        axis.values = seq(0,max(mtx.up.sorted),length.out = 10)
        if(input$cmbMatrix == "matrix_PVAL.txt"){
          mtext('P-value',side=2,line=4)
          axis(side = 2, at = axis.values,labels=scientific_format(1)(1/10^abs(axis.values)),las=1)
        } else{
          mtext('Odds-ratio',side=2,line=5)
          axis(side = 2, at = axis.values, labels=scientific_format(2)(2^axis.values),las=1)
        }
        # draw rotated x-labels
        axis(side=1, labels = FALSE,tick=F)
        text(x.loc, par("usr")[3], srt = 45,  adj=c(1,1),
             labels = head(rownames(mtx.up.sorted),30), xpd = TRUE) 
      }
      # print Enrichmentdown bar plot to PDF
      mtx.down.sorted <- abs(getEnrichmentDown())
      if (nrow(mtx.down.sorted)==0){
        # plot raw text
        par(mar = c(0,0,0,0))
        plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
        text(x = 0.5, y = 0.5, paste("Nothing overrepresented to plot."),  cex = 1.6, col = "black")
        box()
      } else{
        par(mar = c(10,7,4.1,2.1))
        if (nrow(mtx.down.sorted) < 3){
         x.loc <- barplot(as.matrix(t(head(mtx.down.sorted,30))), beside=T,col = "green4",
                  space=c(0.2,1), cex.names=1, las=2, names.arg=head(rownames(mtx.down.sorted),30),xaxt='n',main = "Depleted epigenomic associations",xlim=c(0,10*nrow(mtx.down.sorted)),axes = F)
        }else{
         x.loc <- barplot(as.matrix(t(head(mtx.down.sorted,30))), beside=T,col = "green4",
                  space=c(0.2,1), cex.names=1, las=2, names.arg=head(rownames(mtx.down.sorted),30),xaxt='n',main = "Depleted epigenomic associations", axes = F)
        }
        abline(a=0,b=0)
        axis.values = seq(0,max(mtx.down.sorted),length.out = 10)
        if(input$cmbMatrix == "matrix_PVAL.txt"){
          mtext('P-value',side=2,line=4)
          axis(side = 2, at = axis.values,labels=scientific_format(1)(1/10^abs(axis.values)), las=1)
        }
        else{
          mtext('Odds-ratio',side=2,line=5)
          axis(side = 2, at = axis.values, labels=scientific_format(2)(2^axis.values), las=1)
        }
        # draw rotated x-labels
        axis(side=1, labels = FALSE,tick = F)
        text(x.loc, par("usr")[3], srt = 45, adj=c(1,1),
             labels = head(rownames(mtx.down.sorted),30), xpd = TRUE) 
      }
      
      dev.off()
    },
    contentType = 'application/pdf'
  )
  
  
  # Create a different UI depending if there are multiple GF in the results
  output$mainpage <-renderUI({
    # manually load matrix since controls are not loaded yet
    validate(need(try(mtx <- load_gr_data(paste(get.results.dir(), "matrix_PVAL.txt",sep=""))),"Error loading files. Does the data exist?"))
    file.names.annotation <- list.files(paste(get.results.dir(),"annotations/",sep=""))
    
    single.feature = TRUE;
    if (ncol(mtx)>1 & nrow(mtx)>1){single.feature = FALSE}
    if (single.feature == FALSE){
      tabsetPanel(id="tabsMultiple",
                  tabPanel("Enrichment analysis heatmap",
                           downloadButton('downloadEnrichHeatmap',"Download PDF"),
                           d3heatmapOutput("heatmapEnrich", width = "100%", height = "600px"),
                           plotOutput("legendEnrich",width="300px",height="150px")
                  ), 
                  tabPanel("Enrichment analysis barplot",
                           downloadButton('downloadEnrichBarPDF', 'Download PDF'),
                           plotOutput("pltEnrichUp",width="100%",height = "350px"),
                           plotOutput("pltEnrichDown", width="100%", height= "350px")
                  ),
                  tabPanel("Enrichment analysis tables",
                           br(),br(),
                           downloadButton('downloadEnrichTable', 'Download table in tab-separated format'),
                           br(),br(),
                           DT::dataTableOutput("tblEnrichment")),
                  tabPanel("Epigenetic similarity analysis heatmap",
                           fluidPage(
                             fluidRow(
                               column(6,
                                      downloadButton('downloadEpisimHeatmap', 'Download PDF'),
                                      d3heatmapOutput("heatmapEpisim", width = "100%", height = "600px"),
                                      plotOutput("legendEpisim",width="300px",height="150px")
                               ),
                               column(6, 
                                      plotOutput("pltDend",width = "100%", height = "500px")
                               )
                             )
                           )),
                  tabPanel("Epigenetic similarity analysis tables",
                           selectInput("cmbEpigenetics", "Select which epigenetic analysis to show", choices = list("Results not ready yet.")),
                           br(),br(),
                           downloadButton('downloadEpigenetics', 'Download table in tab-separated format'),
                           br(),br(),
                           DT::dataTableOutput("tblEpigenetics")),
                  if (length(file.names.annotation)>0){
                    tabPanel("Annotation Analysis",
                           DT::dataTableOutput("tblAnnotation"))
                  }
      )
    } else{ # this UI is created when only a single GF result is returned
      tabsetPanel(id="tabsSingleGF",
                  tabPanel("Enrichment analysis barplot",
                           downloadButton('downloadEnrichBarPDF', 'Download PDF'),
                           plotOutput("pltEnrichUp",width="100%",height = "350px"),
                           plotOutput("pltEnrichDown", width="100%", height= "350px")
                  ),
                  tabPanel("Enrichment analysis tables",
                           br(),br(),
                           downloadButton('downloadEnrichTable', 'Download table in tab-separated format'),
                           br(),br(),
                           DT::dataTableOutput("tblEnrichment")),
                  if (length(file.names.annotation)>0){
                    tabPanel("Annotation Analysis",
                             DT::dataTableOutput("tblAnnotation"))
                  }
      )
    }
  })
  output$sidebar <- renderUI({
    # manually load matrix since controls are not loaded yet
    validate(need(try(mtx <- load_gr_data(paste(get.results.dir(), "matrix_PVAL.txt",sep=""))),""))
    mtx.col.names <- colnames(data.frame(mtx)) # R converts '-' to '.' in data frames
    file.names.annotation <- list.files(paste(get.results.dir(),"annotations/",sep=""))
    single.feature = TRUE
    if (ncol(mtx)>1 & nrow(mtx)>1){single.feature = FALSE}
    if (single.feature == FALSE){
      sidebarPanel(width = 4,h3("Data Settings"),
                   selectInput("cmbMatrix", label = "Results to visualize", 
                               choices = list("P-values" = "matrix_PVAL.txt", 
                                              "Odds Ratios" = "matrix_OR.txt")),
                   bsTooltip("cmbMatrix", "Select P-value or Odds ratio", placement = "top", trigger = "hover"),
                   conditionalPanel("input.tabsMultiple == 'Enrichment analysis barplot' || input.tabsMultiple == 'Enrichment analysis tables'",
                                    selectInput("cmbFOI", "Select which SNP set to visualize", choices =   mtx.col.names)
                                    ),
                   conditionalPanel("input.cmbMatrix=='matrix_PVAL.txt'",
                                    selectInput("cmbPvalAdjustMethod",label = "P-value multiple testing correction method",
                                                choices = c( "fdr","none","BH","holm", "hochberg", "hommel", "bonferroni","BY"))),
                   conditionalPanel("input.tabsMultiple == 'Annotation Analysis'",
                                    if (length(file.names.annotation)>0){
                                      selectInput("cmbAnnotation", label = "Annotation results to visualize", 
                                                  choices = file.names.annotation)
                                      }),
                   conditionalPanel("input.tabsMultiple == 'Enrichment analysis heatmap' || input.tabsMultiple == 'Epigenetic similarity analysis heatmap'",
                                    hr(),h3("Visualization option"),
                                    selectInput("cmbClustMethod",label = "Clustering method (hclust)", 
                                                choices = list("ward.D",
                                                               "ward.D2",
                                                               "single",
                                                               "complete",
                                                               "average",
                                                               "mcquitty",
                                                               "median",
                                                               "centroid")
                                    )
                   ),
                   conditionalPanel("input.tabsMultiple == 'Epigenetic similarity analysis heatmap'",
                                    selectInput('cmbEpisimCorType',label = "Correlation coefficient type",
                                                choices = list("Pearson's" = "pearson",
                                                               "Spearman's" = "spearman"))
                   ),
                  
                   conditionalPanel("input.tabsMultiple == 'Epigenetic similarity analysis heatmap'",
                                    hr(),h3("Epigenetic similarity"),
                                    sliderInput("sldEpisimNumClust","Number of clusters",min = 2,max=10,value = 2)
                   )
                   
      )
    }else{ # this is for a single column result file
      sidebarPanel(h3("Global Settings"), hr(),
                   selectInput("cmbMatrix", label = "Results to visualize", 
                               choices = list("P-values" = "matrix_PVAL.txt", 
                                              "Odds Ratios" = "matrix_OR.txt")),
                   selectInput("cmbFOI", "Select which SNP set to visualize", choices =  mtx.col.names),
                   conditionalPanel("input.tabsSingleGF == 'Annotation Analysis'",
                                    if (length(file.names.annotation)>0){
                                      selectInput("cmbAnnotation", label = "Annotation results to visualize", 
                                                  choices = file.names.annotation)
                                    }),
                   conditionalPanel("input.cmbMatrix=='matrix_PVAL.txt'",
                                    if(nrow(mtx)>1){
                                      selectInput("cmbPvalAdjustMethod",label = "P-value multiple testing correction method",
                                                  choices = c( "fdr","none","BH","holm", "hochberg", "hommel", "bonferroni","BY"))}
                   )
      )
    }
  })
})

