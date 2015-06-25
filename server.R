suppressMessages(source("utils2.R")) # See the required packages there
suppressMessages(source("episimilarity.R"))
(source("genomerunner_d3heatmap.R"))
library(d3heatmap)
library(dendextendRcpp) # required for extracting the height from the dendrogram
library(tools)
library(colorRamps)

results.dir <- "/home/lukas/db_2.00_06-10-2015/results/test2_single_col/"
genomerunner.mode <- FALSE
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
    mtx <- load_gr_data(paste(get.results.dir(), input$cmbEnrichHeatmap,sep=""))

  })
  
  output$heatmapEnrich <- renderD3heatmap({
    mat <- get.matrix()
    coloring<-colorRampPalette(c("blue", "yellow", "red"))
  
    if(!check.single_gf()){
      d3heatmap::d3heatmap(as.matrix(mat),heatmap_options = list(hclust=function(tmp) {hclust(tmp, method = input$cmbEnrichClust)}), colors = coloring(coloring.num), show_tip=FALSE,dendro.rds.path=paste(get.results.dir(),"heatmap.dend.rds", sep=""))
    }else{
      NULL
    }
  })
  
  output$legendEnrich <- renderPlot({
    mtx <- get.matrix()
    coloring<-colorRampPalette(c("blue", "yellow", "red"))
    color.range = seq(min(mtx),max(mtx), (max(mtx)-min(mtx))/coloring.num )
    plot(color.range,rep(1,coloring.num+1),col=coloring(coloring.num+1),pch=15,cex=10,main="Heatmap Legend",ylab="",xlab="")
  })
  
  output$tblEnrichment <-renderDataTable({      
    
    enrichment.data <- read.csv(paste(get.results.dir(),"enrichment/",input$cmbEnrichTable,".txt",sep = ""),sep="\t")
    #convert values to numeric form for sorting purposes
    for(x in list(2,4)){
      enrichment.data[[x]] <- as.numeric(enrichment.data[[x]])
    }
    enrichment.data
  })
  
  ## enrichment up and down plots for single column
  get.bar.plot.data <- reactive({
    mtx <- get.matrix() # Make data frame, to allow row names
    ##rownames(mtx) <- mtx$GF; mtx <- mtx[, -1] # Make row names
    # Summarize the values by cell type and factor
    pmax <- function(x) { x[order(abs(x), decreasing=T)][1] } # Get absolute maximum p-value, keeping sign
    ##mtx <- mtx %>% dplyr::select(-description) %>% group_by(cell, factor) 
    mtx <- melt(mtx, id.vars=c("cell", "factor"), , variable.name = "experiment", value.name = "pvalue")
    mtx <- mtx %>% group_by(cell, factor, experiment) %>% dplyr::summarise(pmax = pmax(pvalue))
    mtx <- dcast(mtx, cell + factor ~ experiment, value.var="pmax")
    mtx <- cbind(dplyr::select(mtx, 3:ncol(mtx)), dplyr::select(mtx, cell, factor))
  })
  
  check.single_gf <- reactive({
    # Check if there is only one column in the matrix.  If so, we will plot a bar plots intead of heatmap
    mtx <-get.matrix()
    if(ncol(mtx)==1){
      def.value = 30
      if (nrow(mtx)<30){def.value = nrow(mtx)}
      updateSliderInput(session,"sldNumFeatures",min = 1,max = nrow(mtx),value = def.value)
      return(TRUE)
    }else{
      return(FALSE)
    }
  })
  
  
  output$pltEnrichUp <- renderPlot({
    if (check.single_gf() != TRUE){
      return(plot.new())
    }
    updown.split = switch (input$cmbEnrichHeatmap,"matrix_PVAL.txt" = 0, "matrix_OR.txt" = 1) 
    mtx <-  data.frame(get.matrix())
    mtx.up <- subset(mtx,mtx[1] > updown.split)
    # filter out results that do not meet pvalue threshold
    log10.pval = -log10(input$numBarplotThreshold)
    if (input$cmbEnrichHeatmap == "matrix_PVAL.txt"){
      mtx.up <- subset(mtx.up, mtx.up[1] > log10.pval, drop=F)
      }
    if (nrow(mtx.up)==0){
      # plot raw text
      par(mar = c(0,0,0,0))
      plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
      text(x = 0.5, y = 0.5, paste("Nothing overrepresented to plot."),  cex = 1.6, col = "black")
      box()
      return()
    }
    # sort the results
    mtx.up.sorted <- mtx.up[order(mtx.up[1],decreasing = T),,drop=FALSE]
    
    barplot(as.matrix(t(head(mtx.up.sorted,input$sldNumFeatures))), beside=T,col = "red3",
            space=c(0.2,1), cex.names=0.8, las=2, names.arg=head(rownames(mtx.up.sorted),input$sldNumFeatures),ylab="-log10(p-value)",main="Enriched epigenomic associations")
    abline(a=0,b=0)
    
    #barplot1(head(mtx.up.sorted,input$sldNumFeatures),names.args = head(rownames(mtx.up.sorted),input$sldNumFeatures))
    #names.args.up[names.args.up == "NA:NA"] <- make.names(unlist(lapply(mtx.sorted.up, rownames)), unique=T)[names.args.up == "NA:NA"]
    #names.args.dn[names.args.dn == "NA:NA"] <- make.names(unlist(lapply(mtx.sorted.dn, rownames)), unique=T)[names.args.dn == "NA:NA"]
    #bottom <- 8
    # Plot barplots
    #if (!grepl("dn", toPlot) & (nrow(mtx.barplot.up) > 0)) { barplot1(mtx.barplot.up[, seq(1:length(colnum)), drop=F], "topright", bottom=bottom, names.args=names.args.up, pval=pval) }
  })
  
  output$pltEnrichDown <- renderPlot({
    if (check.single_gf() != TRUE){
      return(plot.new())
    }
    updown.split = switch (input$cmbEnrichHeatmap,"matrix_PVAL.txt" = 0, "matrix_OR.txt" = 1) 
    mtx <-  data.frame(get.matrix())
    mtx.down <- subset(mtx,mtx[1] < updown.split)
    # filter out results that do not meet pvalue threshold
    log10.pval = -log10(input$numBarplotThreshold)
    if (input$cmbEnrichHeatmap == "matrix_PVAL.txt"){
      mtx.down <- subset(mtx.down, mtx.down[1] < -log10.pval, drop=F)
    }
    if (nrow(mtx.down)==0){
      # plot raw text
      par(mar = c(0,0,0,0))
      plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
      text(x = 0.5, y = 0.5, paste("Nothing underrepresented to plot."),  cex = 1.6, col = "black")
      box()
      return()
    }
    mtx.down.sorted <- mtx.down[order(mtx.down[1],decreasing = F),,drop=FALSE]
    barplot(as.matrix(t(head(mtx.down.sorted,input$sldNumFeatures))), beside=T,col = "green4",
            space=c(0.2,1), cex.names=0.8, las=2, names.arg=head(rownames(mtx.down.sorted),input$sldNumFeatures),ylab="-log10(p-value)\nnegative = underrepresentation",main = "Depleted epigenomic associations")
    abline(a=0,b=0)
    #barplot(head(mtx.down.sorted,input$sldNumFeatures),names.args = head(rownames(mtx.down.sorted),input$sldNumFeatures))
  })
  
  # episimilarity ---------------------------------------------------------------
  get.corr.matrix <- reactive({
    mtx <- load_gr_data(paste(get.results.dir(), input$cmbEpisimHeatmap,sep=""))
    mtx <- scale(mtx)
    
    rcorr(as.matrix(mtx), type=input$cmbEpisimCorType)[[1]]
  })
  
  get.cor.hclust.dendrogram <- reactive({
    cor.mat <- get.corr.matrix() 
    as.dendrogram(hclust(as.dist(1-cor.mat), method=input$cmbEpisimClustMethod))
  })
  
  output$heatmapEpisim <- renderD3heatmap({     
    cor.mat <- get.corr.matrix()
    hclustergram <- get.cor.hclust.dendrogram()
    coloring<-colorRampPalette(c("blue", "yellow", "red"))
    d3heatmap::d3heatmap(as.matrix(cor.mat),heatmap_options = list(Rowv=hclustergram,Colv=hclustergram,keep.dendro=TRUE),colors = coloring(coloring.num),dendro.rds.path=paste(get.results.dir(),"heatmap.dend.rds", sep=""))
  })
  
  output$legendEpisim <- renderPlot({
    mtx <- get.corr.matrix()
    coloring<-colorRampPalette(c("blue", "yellow", "red"))
    color.range = seq(min(mtx),max(mtx), (max(mtx)-min(mtx))/coloring.num )
    plot(color.range,rep(1,coloring.num+1),col=coloring(coloring.num+1),pch=15,cex=10,main="Heatmap Legend",ylab="",xlab="")
  })
  
  output$tblEpigenetics <-renderDataTable({
    mtx.deg <- readRDS(file=paste(get.results.dir(),"mtx.deg.episim.RDS",sep = "/"))
    
    #convert values to numeric form for sorting purposes
    for(x in list("adj.p.val",3,4)){
      mtx.deg[[input$cmbEpisimTable]][[x]] <- as.numeric(mtx.deg[[input$cmbEpisimTable]][[x]])
    }
    mtx.deg[[input$cmbEpisimTable]]
  })
  
  output$pltDend <- renderPlot({ 
    cor.mat <- get.corr.matrix() # this line ensure that dendrogram is redrawn when heatmap is
    hclustergram <- get.cor.hclust.dendrogram() # ensures that dendrogram is redrawn when hclust method is changed
    
    dend = readRDS(file = paste(get.results.dir(), "heatmap.dend.rds",sep=""))
    plot(as.dendrogram(dend, hang=-1)) # Plot dendrogram
    cl_num <- input$sldEpisimNumClust # Empirically set desired numter of clusters
    cols <- rainbow(cl_num) # Make colos for clusters
    print(cl_num)
    hcut <- heights_per_k.dendrogram(dend)[cl_num] # extract the height to cut based on # of groups
    # get the cluster labels
    mtx.clust <- dend %>% mtx.clusters(height=hcut, minmembers=3)
    # write.table(as.data.frame(mtx.clust), "/home/lukas/clustering_all.txt", sep="\t", row.names=FALSE, quote=FALSE)
    mtx = load_gr_data(paste(get.results.dir(), input$cmbEpisimHeatmap,sep="")) # load the original matrix
    mtx.deg <- mtx.degfs(mtx[, mtx.clust$eset.labels], mtx.clust, label="broadPeak2")
    
    updateSelectInput(session,"cmbEpigenetics","Select which epigenetic table to render",choices = names(mtx.deg))
    mtx.deg.path = paste(get.results.dir(),"mtx.deg.episim.RDS",sep = "/")
    saveRDS(mtx.deg,file=mtx.deg.path)
    # create the tabs
    
    rect.hclust(as.hclust(dend), k=cl_num, border=cols) # Define the clusters by rectangles
    
  })
})
