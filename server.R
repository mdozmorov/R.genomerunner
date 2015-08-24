# suppressMessages(source("utils2.R")) # See the required packages there
# suppressMessages(source("episimilarity.R"))
source("functions/load_required_packages.R")
source("functions/load_gr_data.R")
source("functions/mtx.transform.R")
source("functions/mtx.untransform.R")
source("functions/mtx.clusters.R")
source("functions/mtx.degfs.R")
source("functions/mtx.cellspecific.R")
#shiny::runApp(host='0.0.0.0',port=4494)

results.dir <- "/home/lukas/db_2.00_06-10-2015/results/largerun/"
# Mikhail paths
# results.dir <- "/home/mdozmorov/db_5.00_07-22-2015/results/"
results.dir <- "/Users/mikhail/Documents/tmp/results/py4rdowe2jhj90jd9y1sy4c0gsz7t9bj/"

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
    mtx
  })
  
  output$heatmapEnrich <- renderD3heatmap({
    withProgress({
      # force heatmap to be redrawn when controls change
      untransform.method <- "none"
      # mtx <- get.adjust.matrix()
      if (input$cmbMatrix == "matrix_PVAL.txt"){
        mtx <- get.adjust.matrix()
        untransform.method <- "log10"
      }else{
        mtx <- get.matrix()
        untransform.method <- 'log2'
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
      dend.path <- paste(get.results.dir(),"heatmap.dend.rds", sep="")
      par(cex.main=0.65, oma=c(2,0,0,5), mar=c(5, 4.1, 4.1, 5)) # Adjust margins
      coloring<-colorRampPalette(c("blue", "yellow", "red"))
     
      d3heatmap::d3heatmap(as.matrix(mtx),hclust=function(tmp) {hclust(tmp, method = input$cmbClustMethod)}, colors = coloring(coloring.num), tip_transformation = untransform.method,
                         xaxis_font_size = "10pt", yaxis_font_size = "10pt",xaxis_height=200,yaxis_height=200,dendro.rds.path=dend.path)
      }, message = "Rendering heatmap",value = 1.0)
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
      
      mtx.table <- cbind(scientific_format(3)(abs(mtx[selectedFOI])), pval.sig, scientific_format(3)(mtx.adjust))
      colnames(mtx.table) <- c("p.value","direction","adj.p.val")
    }else{ 
      # for odds ratio
      or.sig <- rep("Not significant", nrow(mtx))
      or.sig[mtx[, selectedFOI] > 1] <- "Enriched"
      or.sig[mtx[, selectedFOI] < 1] <- "Depleted"
      mtx.table <- cbind(abs(mtx[selectedFOI]),
                         or.sig)
      colnames(mtx.table) <- c("odds.ratio","direction")
    }
    mtx.table = cbind(GF.Name = rownames(mtx.table),mtx.table) # make the GF.name a column instead of just the rowname so we can left_join
    rownames(mtx.table) <- NULL
    # gfAnnot is loaded in utils2
    
    mtx.table <- left_join(mtx.table,gfAnnot[,(names(gfAnnot) %in% c("file_name", "cell", "cell_desc", "factor", "factor_desc", "source", "source_desc"))],
                           by=c("GF.Name"="file_name"))
    
    return(mtx.table)
  })
  
  output$tblEnrichment <-renderDataTable({
    withProgress({
      num.char <- 50
      table.enrich <- get.enrichment.table()
      table.enrich <- apply(table.enrich,c(1,2),function(x) {
        if (!is.na(x) & nchar(x)>num.char){
          return(paste(substring(x,1,num.char),  "..."))
        } else{
          return(x)
        }
      })
    }, message = "Loading enrichment table",value=1.0)
  }, options = list( lengthMenu = list(c(10, 50, 100,-1), c('10', '50','100', 'All')),
                     pageLength = 10))
  
  
  
  output$downloadEnrichTable <- downloadHandler(
    filename = function() { 
      return("Enrichment_table.txt")
    },
    content = function(file) {
      write.table(x = get.enrichment.table(),file =  file ,sep = "\t",quote = F,row.names = F)
    }
  )
  
  get.annotation.table <- reactive({
    withProgress({
    validate(need(try(mtx <- read.table(paste(get.results.dir(),"annotations/", input$cmbAnnotation,sep=""),header=T)),"Annotation results not available."))
    mtx}, message = "Loading Annotation table",value=1.0)
  })
  
  output$tblAnnotation <- renderDataTable({
    get.annotation.table()
  },options = list( lengthMenu = list(c(10, 50, 100,-1), c('10', '50','100', 'All')),
                     pageLength = 10))
  
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
              space=c(0.2,1), cex.names=1, las=2, xaxt='n',main="Enriched regulatory associations",xlim=c(0,10*nrow(mtx.up.sorted)),axes=F)
    } else{
     x.loc =  barplot(as.matrix(t(head(mtx.up.sorted,30))), beside=T,col = "red3",
              space=c(0.2,1), cex.names=1, las=2, xaxt='n',main="Enriched regulatory associations",axes=F)
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
              space=c(0.2,1), cex.names=1, las=2, names.arg=head(rownames(mtx.down.sorted),30),xaxt='n',main = "Depleted regulatory associations",xlim=c(0,10*nrow(mtx.down.sorted)),axes = F)
    }else{
      x.loc = barplot(as.matrix(t(head(mtx.down.sorted,30))), beside=T,col = "green4",
              space=c(0.2,1), cex.names=1, las=2, names.arg=head(rownames(mtx.down.sorted),30),xaxt='n',main = "Depleted regulatory associations", axes = F)
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
    # If there are columns with SD=0, add jitter to it. Needed for pair-wise column correlation analysis (regulatory similarity analysis). Only valid if there's more than 1 row
    if (nrow(mtx) > 1) {
      ind <- apply(mtx, 2, function(x) sd(x, na.rm=TRUE)) == 0 # Indexes of such columns
      if (sum(ind) > 0) {
        set.seed(1)
        mtx[, ind] <- jitter(mtx[, ind, drop=FALSE], factor=0.1)
      }
    }
    mtx <- scale(mtx)
    
   
   mtx <- rcorr(as.matrix(mtx), type=input$cmbEpisimCorType)[[1]]
   write.table(x = mtx,file = paste(results.dir,"matrix_CORR.txt",sep=""))
   mtx
  })
  
  get.cor.hclust.dendrogram <- reactive({
    cor.mat <- get.corr.matrix() 
    as.dendrogram(hclust(as.dist(1-cor.mat), method=input$cmbClustMethod))
  })
  
  output$heatmapEpisim <- renderD3heatmap({ 
    withProgress({
      mat <- get.matrix()
      validate(need(ncol(mat)>2,"Need at least 3 SNPs of interest files to perform clustering."))
      validate(need(nrow(mat)>4,"Need at least 5 genome features to perform clustering."))
      cor.mat <- get.corr.matrix()
      hclustergram <- get.cor.hclust.dendrogram()
      coloring<-colorRampPalette(c("blue", "yellow", "red"))
      # TODO: Colv == Rowv?
      d3heatmap::d3heatmap(as.matrix(cor.mat),Rowv=hclustergram,Colv='Rowv',colors = coloring(coloring.num),dendro.rds.path=paste(get.results.dir(),"heatmap.dend.rds", sep=""),
                           xaxis_font_size = "10pt", yaxis_font_size = "10pt",xaxis_height=200,yaxis_height=200)
    }, message = "Rendering heatmap",value = 1.0)
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
    mtx.deg <- suppressWarnings(mtx.degfs(mtx[, mtx.clust$eset.labels], mtx.clust, isOR = is.OR))
    
    updateSelectInput(session,"cmbEpigenetics","Select which comparison to show",choices = names(mtx.deg))
    mtx.deg.path = paste(get.results.dir(),"mtx.deg.episim.RDS",sep = "")
    saveRDS(mtx.deg,file=mtx.deg.path)
    return(mtx.deg)
  })
  
  get.epigenetics.table <- reactive({
    withProgress({
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
    }, message = "Calculating regulatory differences",value = 1.0)
    #convert values to numeric form for sorting purposes
    if(input$cmbMatrix == "matrix_PVAL.txt"){
      for(x in list("adj.p.val",3,4)){
        mtx.deg[[selectedCor]][[x]] <- scientific_format(3)(as.numeric(mtx.deg[[selectedCor]][[x]]))
      }
    } else {
        mtx.deg[[selectedCor]][["adj.p.val"]] <- scientific_format(3)(as.numeric(mtx.deg[[selectedCor]][["adj.p.val"]]))
      }
    mtx.deg[[selectedCor]][, !(colnames(mtx.deg[[1]]) %in% c("full_path", "URL", "full_description", "category", "category_desc"))]
  })
  
  output$tblEpigenetics <-renderDataTable({
    mtx <- get.matrix()
    validate(need(ncol(mtx)>2,"Need at least 3 SNPs of interest files to perform clustering."))
    validate(need(nrow(mtx)>4,"Need at 5 least genome features to perform clustering."))
    validate(need(try(get.epigenetics.table()),"Either nothing is significant, or there are too few SNP sets per cluster. Re-run the analysis using more SNP sets, or try a different clustering method."))
    table.epi <- get.epigenetics.table()
    num.char <- 50
    table.epi <- apply( table.epi,c(1,2),function(x) {
      if (!is.na(x) & nchar(x)>num.char){
        return(paste(substring(x,1,num.char),  "..."))
      } else{
        return(x)
      }
    })
  },options = list( lengthMenu = list(c(10, 50, 100,-1), c('10', '50','100', 'All')),
                    pageLength = 10))
  
  output$downloadEpigenetics <- downloadHandler(
    filename = function() { 
      return("Regulatory_table.txt")
    },
    content = function(file) {
      write.table(x = get.epigenetics.table(),file =  file ,sep = "\t",quote = F,row.names = F)
    }
  )
  
  output$pltDend <- renderPlot({ 
    withProgress({
      mtx <- get.matrix()
      validate(need(ncol(mtx)>2,""))
      validate(need(nrow(mtx)>4,""))
      
      cor.mat <- get.corr.matrix() # this line ensure that dendrogram is redrawn when heatmap is
      hclustergram <- get.cor.hclust.dendrogram() # ensures that dendrogram is redrawn when hclust method is changed
      
      dend = readRDS(file = paste(get.results.dir(), "heatmap.dend.rds",sep=""))
      plot(as.dendrogram(dend, hang=-1)) # Plot dendrogram
      cl_num <- input$sldEpisimNumClust # Empirically set desired numter of clusters
      cols <- rainbow(cl_num) # Make colos for clusters
      hcut <- heights_per_k.dendrogram(dend)[cl_num] # extract the height to cut based on # of groups
      # get the cluster labels
      mtx.clust <-validate(need(try(dend %>% mtx.clusters(height=hcut, minmembers=3)),"Try using a lower number of clusters"))
      # write.table(as.data.frame(mtx.clust), "/home/lukas/clustering_all.txt", sep="\t", row.names=FALSE, quote=FALSE)
      mtx = load_gr_data(paste(get.results.dir(), input$cmbMatrix,sep="")) # load the original matrix
      
      # Define the clusters by rectangles
      validate(need(try(rect.hclust( as.hclust(dend), k=cl_num, border=cols)),"Try a different clustering method."))
    }, message = "Rendering dendrogram", value = 1.0)
      
  })
  
  calculateCTEnrichment <- reactive({
    mtx <- load_gr_data(paste(get.results.dir(), 'matrix_PVAL.txt',sep=""))
    validate(need(nrow(mtx)>5,"Insufficient data for performing cell type-specific enrichment analysis"))
    #running function
    mtx.CTE <- mtx.cellspecific(mtx)
  })
  
  get.CTEnrichment.table <- reactive({
    mtx.CTE <- calculateCTEnrichment()
    if (is.character(mtx.CTE)) {
      return(data.frame(NoResults="Insufficient data for performing cell type-specific enrichment analysis"))
    }
    selectedCor <- input$cmbFOI
    if (is.character(mtx.CTE[[selectedCor]])) {
      return(data.frame(NoResults="Nothing significant"))
    }
    return(mtx.CTE[[selectedCor]])
  })
  
  output$tblCTEnrichment <- renderDataTable({
    withProgress({
      table.CTE <-  get.CTEnrichment.table()
      num.char <- 50
      table.CTE  <- apply( table.CTE ,c(1,2),function(x) {
        if (!is.na(x) & nchar(x)>num.char){
          return(paste(substring(x,1,num.char),  "..."))
        } else{
          return(x)
        }
      })
    }, message = "Loading table", value = 1.0)
  },options = list( lengthMenu = list(c(10, 50, 100,-1), c('10', '50','100', 'All')),
                     pageLength = 10))
  
  output$downloadCTEnrichment <- downloadHandler(
    filename = function() { 
    return("EnrichmentCT_table.txt")
  },
  content = function(file) {
    write.table(x = get.CTEnrichment.table(),file =  file ,sep = "\t",quote = F,row.names = T,col.names=NA)
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
      par(cex.main=0.65, oma=c(2,0,0,5), mar=c(5, 4.1, 4.1, 5)) # Adjust margins
      heatmap.2(as.matrix(mtx),hclust=function(tmp) {hclust(tmp, method = input$cmbClustMethod)},
                trace="none", density.info="none", col=colorRampPalette(c("blue", "yellow", "red")), main = "Enrichment Heatmap",
                cexRow=0.8, cexCol=1, margins = c(10,10), srtRow = 0, srtCol = 45)
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
      par(cex.main=0.65, oma=c(2,0,0,5), mar=c(5, 4.1, 4.1, 5)) # Adjust margins
      coloring<-colorRampPalette(c("blue", "yellow", "red"))
      heatmap.2(as.matrix(cor.mat), Colv = hclustergram, Rowv = hclustergram,
                trace="none", density.info="none", col=colorRampPalette(c("blue", "yellow", "red")), main = "Episimilarity Heatmap",
                cexRow=0.8, cexCol=1, margins = c(10,10), srtRow = 0, srtCol = 45)
      dev.off()
    },
    contentType = 'application/pdf'
  )
  
  output$downloadAnnotation <- downloadHandler(
    filename = function() { 
      return("Enrichment_table.txt")
    },
    content = function(file) {
      write.table(x = get.annotation.table(),file =  file ,sep = "\t",quote = F,row.names = F)
    }
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
                  space=c(0.2,1), cex.names=1, las=2, names.arg=head(rownames(mtx.up.sorted),30),xaxt='n',main="Enriched regulatory associations",xlim=c(0,10*nrow(mtx.up.sorted)),axes=F)
        } else{
          x.loc <-barplot(as.matrix(t(head(mtx.up.sorted,30))), beside=T,col = "red3",
                  space=c(0.2,1), cex.names=1, las=2, names.arg=head(rownames(mtx.up.sorted),30),xaxt='n',main="Enriched regulatory associations",axes=F)
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
                  space=c(0.2,1), cex.names=1, las=2, names.arg=head(rownames(mtx.down.sorted),30),xaxt='n',main = "Depleted regulatory associations",xlim=c(0,10*nrow(mtx.down.sorted)),axes = F)
        }else{
         x.loc <- barplot(as.matrix(t(head(mtx.down.sorted,30))), beside=T,col = "green4",
                  space=c(0.2,1), cex.names=1, las=2, names.arg=head(rownames(mtx.down.sorted),30),xaxt='n',main = "Depleted regulatory associations", axes = F)
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
  
  output$downloadZIP <- downloadHandler(
    filename = function() { 
      return("genome_runner.zip")
    },
    content = function(file) {
      # ensure correlation matrix is created
      mtx <- get.matrix()
      if (ncol(mtx)>1 & nrow(mtx)>1){get.corr.matrix()}
      
      
      # Append gfAnnot columns to the end of the PVAL and OR matrix
      mtx <- read.csv(paste(get.results.dir(),"matrix_PVAL.txt",sep=""),sep = '\t')
      mtx <- data.frame(GF=rownames(mtx), mtx)
      mtx <- left_join(mtx, gfAnnot, by=c("GF" = "file_name"))
      rownames(mtx) <- mtx$GF; mtx$GF <- NULL
      write.table(mtx,file=paste(get.results.dir(),"matrix_PVAL_annot.txt",sep=""), sep="\t", quote=FALSE, col.names=NA)
      
      mtx <- read.csv(paste(get.results.dir(),"matrix_OR.txt",sep=""),sep = '\t')
      mtx <- data.frame(GF=rownames(mtx), mtx)
      mtx <- left_join(mtx, gfAnnot, by=c("GF" = "file_name"))
      rownames(mtx) <- mtx$GF; mtx$GF <- NULL
      write.table(mtx,file=paste(get.results.dir(),"matrix_OR_annot.txt",sep=""), sep="\t", quote=FALSE, col.names=NA)
      mtx.clust = tryCatch({
         calculate.clust()
       }, error = function(e) {
        "Try different clustering settings"
       })
      
      # save each new data frame as an individual .csv file based on its name
      
      lapply(1:length(mtx.clust), function(i) write.table(mtx.clust[[i]], 
                                                      file = paste0(get.results.dir(),names(mtx.clust[i]), ".txt"),
                                                      quote=F,row.names = F,sep='\t'))
     
      # zip up text files
      files.txt <- as.character(sapply(c("gr_log.txt","detailed.txt","matrix_PVAL_annot.txt","matrix_OR_annot.txt","matrix_CORR.txt"), function(x){
        paste(get.results.dir(),x,sep = "")
      }))
      files.txt <- append(files.txt, paste0(get.results.dir(),names(mtx.clust), ".txt"))
      file.names.annotation <- list.files(paste(get.results.dir(),"annotations",sep=""),full.names = T)
      anot.path <- paste(get.results.dir(),"annotations.zip",sep="")
      if (length(file.names.annotation) != 0 & !file.exists(anot.path)){
        zip(zipfile <- anot.path ,files =  file.names.annotation,flags = "-j")
      }
      
      enrich.path <- paste(get.results.dir(),"enrichment.zip",sep="")
      file.names.enrichment <- list.files(paste(get.results.dir(),"enrichment",sep=""),full.names = T)
      if (length(file.names.enrichment) != 0 & !file.exists(enrich.path)){
        zip(zipfile <- enrich.path,files = file.names.enrichment,flags = "-j")
      }
      file.all <- c(enrich.path,anot.path,files.txt)
      file.all <- gsub("//", "/" ,file.all)
      zip(zipfile <- file,files = file.all,flags = "-j")
    },
    contentType = "application/zip"
  )
  
  
  # Create a different UI depending if there are multiple GF in the results
  output$mainpage <-renderUI({
    # manually load matrix since controls are not loaded yet
    validate(need(try(mtx <- load_gr_data(paste(get.results.dir(), "matrix_PVAL.txt",sep=""))),"Error loading files. Either no significant results are available, or data files are corrupted. Please, re-run the analysis using larger number of genome annotation datasets."))
    file.names.annotation <- list.files(paste(get.results.dir(),"annotations/",sep=""))
    
    single.feature = TRUE;
    if (ncol(mtx)>1 & nrow(mtx)>1){single.feature = FALSE}
    if (single.feature == FALSE){
      tabsetPanel(id="tabsMultiple",
                  tabPanel("Enrichment heatmap",
                           br("Heatmap of the enrichment analysis results. Rows - names of regulatory elements, columns - names of SNP sets, cells - enrichment p-values/odds ratios. Blue/red gradient highlights depleted/enriched associations for corresponding regulatory elements and SNP sets, respectively"),
                           br("Mouse over the heatmap to see numerical values. Click-and-drag to zoom in, single click to reset zoom."),
                           downloadButton('downloadEnrichHeatmap',"Download PDF"),
                           d3heatmapOutput("heatmapEnrich", width = "100%", height = "600px"),
                           plotOutput("legendEnrich",width="300px",height="150px")
                  ), 
                  tabPanel("Enrichment barplot",
                           br("SNP set-specific enrichment results. Height of bars corresponds to the strength of enriched/depleted associations, top/bottom barplot, respectively."),
                           br(),
                           downloadButton('downloadEnrichBarPDF', 'Download PDF'),
                           plotOutput("pltEnrichUp",width="100%",height = "350px"),
                           plotOutput("pltEnrichDown", width="100%", height= "350px")
                  ),
                  tabPanel("Enrichment tables",
                           br("Enrichment analysis results in text format."),
                           br(),
                           downloadButton('downloadEnrichTable', 'Download table'),
                           br(),br(),
                           DT::dataTableOutput("tblEnrichment")),
                  tabPanel("Regulatory similarity heatmap",
                           fluidPage(
                             fluidRow(
                               column(6,
                                      br("Heatmap of regulatory similarity among SNP sets. Cells show correlation coefficients for each pair-wise correlation of SNP set-specific regulatory profiles."),
                                      br("Mouse over the heatmap to see numerical values. Click-and-drag to zoom in, single click to reset zoom."),
                                      downloadButton('downloadEpisimHeatmap', 'Download PDF'),
                                      d3heatmapOutput("heatmapEpisim", width = "100%", height = "600px"),
                                      plotOutput("legendEpisim",width="300px",height="150px")
                               ),
                               column(6,
                                      br("Dedrogram of regulatory similarity among SNP sets. Clusters of SNP sets having strong regulatory similarity indicates these SNP sets are enriched in similar regulatory elements."),
                                      br("Adjust the number of clusters to identify regulatory differences among them on the \"Differential regulatory analysis\" tab."),
                                      plotOutput("pltDend",width = "100%", height = "500px")
                               )
                             )
                           )),
                  tabPanel("Differential regulatory analysis",
                           br("Differential regulatory analysis identifies regulatory elements differentially enriched between clusters of SNP sets (e.g., \"cX_vs_cY\"). Adjust the number of clusters and other clustering metrics on the \"Regulatory similarity heatmap\" tab."),
                           br(),
                           selectInput("cmbEpigenetics", "Select which comparison to show", choices = list("Results not ready yet.")),
                           br(),
                           downloadButton('downloadEpigenetics', 'Download table'),
                           br(),br(),
                           DT::dataTableOutput("tblEpigenetics")),
                  if (length(file.names.annotation)>0){
                    tabPanel("Annotation Analysis",
                           br("Annotation analysis tables. Each SNP in a set (rows) is annotated for overlap with regulatory elements (columns). A non-zero value indicates a SNP overlaps corresponding regulatory element."),
                           br("If more than 100 regulatory elements were selected, the annotation tables are split into multiple tables, each having 100 columns or less."),
                           br(),
                           downloadButton('downloadAnnotation', 'Download Table'),
                           br(),br(),
                           DT::dataTableOutput("tblAnnotation"))
                  }else{
                    conditionalPanel('False',tabPanel("Annotation Analysis")
                    )
                  },
                  tabPanel("Cell-type enrichment analysis",
                           br("Cell-type enrichment analysis detects cell type specificity of the enrichments of SNP sets. It tests whether enrichments in cell type-specific regulatory elements (AvPvalCell) are significantly different from the average enrichments (AvPvalTot)."),
                           br("This analysis requires several enrichment analyses per cell type. E.g., selecting \"DNase\" regulatory elements, with one regulatory set per cell type will provide insufficient information for cell type-specific enrichment analysis. Select categories of regulatory elements having multiple cell type-specific regulatory data, e.g., \"Histone\" and/or \"chromStates\"."),
                           br(),
                           downloadButton('downloadCTEnrichment', 'Download table'),
                           br(),br(),
                           DT::dataTableOutput("tblCTEnrichment")),
                  tabPanel("Download",
                           br(),
                           h3("Download results"),
                           downloadButton('downloadZIP', 'Download enrichment and/or annotation analysis results'))
      )
    } else{ # this UI is created when only a single GF result is returned
      tabsetPanel(id="tabsSingleGF",
                  tabPanel("Enrichment barplot",
                           br("SNP set-specific enrichment results. Height of bars corresponds to the strength of enriched/depleted associations, top/bottom barplot, respectively."),
                           br(),
                           downloadButton('downloadEnrichBarPDF', 'Download PDF'),
                           plotOutput("pltEnrichUp",width="100%",height = "350px"),
                           plotOutput("pltEnrichDown", width="100%", height= "350px")
                  ),
                  tabPanel("Enrichment tables",
                           br("Enrichment analysis results in text format."),
                           br(),
                           downloadButton('downloadEnrichTable', 'Download table'),
                           br(),br(),
                           DT::dataTableOutput("tblEnrichment")),
                  if (length(file.names.annotation)>0){
                    tabPanel("Annotation Analysis",
                             br("Annotation analysis tables. Each SNP in a set (rows) is annotated for overlap with regulatory elements (columns). A non-zero value indicates a SNP overlaps corresponding regulatory element."),
                             br("If more than 100 regulatory elements were selected, the annotation tables are split into multiple tables, each having 100 columns or less."),
                             br(),
                             downloadButton('downloadAnnotation', 'Download Table'),
                             DT::dataTableOutput("tblAnnotation"))
                  }else{
                    conditionalPanel('False',tabPanel("Annotation Analysis")
                    )
                  },
                  tabPanel("Cell-type enrichment analysis",
                           br("Cell-type enrichment analysis detects cell type specificity of the enrichments of SNP sets. It tests whether enrichments in cell type-specific regulatory elements (AvPvalCell) are significantly different from the average enrichments (AvPvalTot)."),
                           br("This analysis requires several enrichment analyses per cell type. E.g., selecting \"DNase\" regulatory elements, with one regulatory set per cell type will provide insufficient information for cell type-specific enrichment analysis. Select categories of regulatory elements having multiple cell type-specific regulatory data, e.g., \"Histone\" and/or \"chromStates\"."),
                           br(),
                           downloadButton('downloadCTEnrichment', 'Download table'),
                           br(),br(),
                           DT::dataTableOutput("tblCTEnrichment")),
                  tabPanel("Download",
                           br(),
                           h3("Download results"),
                           downloadButton('downloadZIP', 'Download enrichment and/or annotation analysis results'))
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
                   conditionalPanel("input.tabsMultiple != 'Cell-type enrichment analysis' && input.tabsMultiple != 'Download'  && input.tabsMultiple != 'Annotation Analysis'",
                     selectInput("cmbMatrix", label = "Results to visualize", 
                                 choices = list("P-values" = "matrix_PVAL.txt", 
                                                "Odds Ratios" = "matrix_OR.txt")),
                     bsTooltip("cmbMatrix", "Select significance or effect size", placement = "right", trigger = "hover")
                   ),
                   conditionalPanel("input.tabsMultiple == 'Enrichment barplot' || input.tabsMultiple == 'Enrichment tables' || input.tabsMultiple == 'Cell-type enrichment analysis'",
                                    selectInput("cmbFOI", "Select which SNP set to visualize", choices =   mtx.col.names)
                                    ),
                   conditionalPanel("input.cmbMatrix=='matrix_PVAL.txt' && input.tabsMultiple != 'Cell-type enrichment analysis' && input.tabsMultiple != 'Differential regulatory analysis'  && input.tabsMultiple != 'Download' && input.tabsMultiple != 'Annotation Analysis'",
                                    selectInput("cmbPvalAdjustMethod",label = "P-value multiple testing correction method",
                                                choices = c( "fdr","none","BH","holm", "hochberg", "hommel", "bonferroni","BY"))),
                   conditionalPanel("input.tabsMultiple == 'Annotation Analysis'",
                                    if (length(file.names.annotation)>0){
                                      selectInput("cmbAnnotation", label = "Annotation results to visualize", 
                                                  choices = file.names.annotation)
                                      }),
                   conditionalPanel("input.tabsMultiple == 'Enrichment heatmap' || input.tabsMultiple == 'Regulatory similarity heatmap'",
                                    hr(),h3("Visualization option"),
                                    selectInput("cmbClustMethod",label = "Clustering method", 
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
                   bsTooltip("cmbClustMethod", "Select clustering method", placement = "right", trigger = "hover"),
                   conditionalPanel("input.tabsMultiple == 'Regulatory similarity heatmap'",
                                    selectInput('cmbEpisimCorType',label = "Correlation coefficient type",
                                                choices = list("Pearson's" = "pearson",
                                                               "Spearman's" = "spearman"))
                   ),
                  
                   conditionalPanel("input.tabsMultiple == 'Regulatory similarity heatmap'",
                                    hr(),h3("Regulatory similarity"),
                                    sliderInput("sldEpisimNumClust","Number of clusters",min = 2,max=10,value = 2)
                   )
      )
    }else{ # this is for a single column result file
      sidebarPanel(h3("Global Settings"), hr(),
                   conditionalPanel("input.tabsSingleGF != 'Cell-type enrichment analysis' && input.tabsSingleGF != 'Download' && input.tabsSingleGF != 'Annotation Analysis'",
                     selectInput("cmbMatrix", label = "Results to visualize", 
                                 choices = list("P-values" = "matrix_PVAL.txt", 
                                                "Odds Ratios" = "matrix_OR.txt"))
                   ),
                   selectInput("cmbFOI", "Select which SNP set to visualize", choices =  mtx.col.names),
                   conditionalPanel("input.tabsSingleGF == 'Annotation Analysis'",
                                    if (length(file.names.annotation)>0){
                                      selectInput("cmbAnnotation", label = "Annotation results to visualize", 
                                                  choices = file.names.annotation)
                                    }),
                   conditionalPanel("input.cmbMatrix=='matrix_PVAL.txt' && input.tabsSingleGF != 'Cell-type enrichment analysis' && input.tabsSingleGF != 'Download' && input.tabsSingleGF != 'Annotation Analysis'",
                                    if(nrow(mtx)>1){
                                      selectInput("cmbPvalAdjustMethod",label = "P-value multiple testing correction method",
                                                  choices = c( "fdr","none","BH","holm", "hochberg", "hommel", "bonferroni","BY"))}
                   )
      )
    }
  })
})

