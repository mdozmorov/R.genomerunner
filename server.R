library(MDmisc)
library(shinyBS)
library(dplyr)
library(ggvis)
library(tidyr)
#shiny::runApp(host='0.0.0.0',port=4494)

#results.dir <- "/home/lukas/db_2.00_06.14.2016/results/"
# Mikhail paths
results.dir <- "/home/lukas/Sample_runs/wo5hpnesw948o05v84lpq05lx6yyaarp/"
#results.dir <- "/home/lukas/db_2.00_06.14.2016/results/gr_ADME_rdmHistone_bPk-processed/"
# results.dir <- "/Users/mikhail/Documents/tmp/results/diseases_vs_rdmHistone_gPk-imputed/"
# results.dir <- "/Users/mikhail/Documents/tmp/results/example2/"
# 
genomerunner.mode <- F
coloring.num = 50
num.char <- 50

shinyServer(function(input, output,session) {

# Parse the GET query string
get.results.dir <- reactive({
  if (genomerunner.mode) {
    query <- parseQueryString(session$clientData$url_search)
    return(paste(results.dir, query$job_id, "/", sep = ""))
  } else {
    return(results.dir)
  }
})
  

# enrichment ------------------------------------------------------------------
get.matrix <- reactive({
  # populate the enrichment table combobox
  file.names.enrichment <- tools::file_path_sans_ext(
    list.files(paste(get.results.dir(), "enrichment/", sep = "")))
  mtx <- gr_load_data(paste(get.results.dir(), input$cmbMatrix, sep = ""))
})
  

get.adjust.matrix <- reactive({
  mtx <- get.matrix()
  sign.mtx <- apply(mtx, 2, sign) # Keep the sign
  adjust.rownames <- rownames(mtx)
  total.columns <- ncol(mtx)
  # Join with annotations
  mtx <- data.frame(GF = adjust.rownames, mtx)  # Attach GF names
  mtx <- dplyr::left_join(mtx, gfAnnot[, c("file_name", "cell")], by = c(GF = "file_name"))
  row.names(mtx) <- adjust.rownames
  class(mtx$cell) <- "character"
  # If some file names is not in the gfAnnot dataframe (e.g., user-provided data), 
  # 'cell' column will contain NAs. replace them with dummy text to allow FDR correction
  mtx$cell[is.na(mtx$cell)] <- "dummy_cell"  
  unique.cells <- unique(mtx$cell)  # Keep unique cell types
  # Adjust for multiple testing and -log10 transform, if needed Process
  # each column Adjust for multiple testing on per-cell-type basis If the
  # cell-specific subset have >1 row, perform correction for multiple
  # testing
  for (i in 1:total.columns) {
    for (u.c in unique.cells) {
      if (sum(mtx$cell == u.c) > 1) {
        # i+1 because we added the GF column
        mtx[mtx$cell == u.c, i + 1] <- 
          gr_transform(apply(gr_untransform(mtx[mtx$cell == u.c, i + 1, drop = F]), 2, 
                             function(x) p.adjust(abs(x), method = input$cmbPvalAdjustMethod)))
      }
    }
  }
  mtx <- sign.mtx * mtx[, 2:(total.columns + 1)] # 'Reattach' the sign
  mtx
})
  

output$heatmapEnrich <- renderD3heatmap({
  withProgress({
    # force heatmap to be redrawn when controls change
    untransform.method <- "none"
    # mtx <- get.adjust.matrix()
    if (input$cmbMatrix == "matrix_PVAL.txt") {
      mtx <- get.adjust.matrix()
      untransform.method <- "log10"
    } else {
      mtx <- get.matrix()
      untransform.method <- "log2"
    }
    n_limit = 30
    # if number of rows > n_limit, sort rows by SD and take top n_limit most variable.
    if (nrow(mtx) > n_limit) {
      mtx.sd <- apply(mtx, 1, sd)
      mtx.sd <- data.frame(mtx.sd)
      # sort rows based on SD
      mtx.sd.order <- mtx[order(mtx.sd, decreasing = T), ]
      mtx.sd.order <- mtx.sd.order[1:n_limit, ]
      mtx <- mtx.sd.order
    }
  }, message = "Rendering heatmap", value = 1)
  # Dendrogram file is created by d3Heatmap
  dend.path <- paste(get.results.dir(), "heatmap.enrich.dend.rds", sep = "")  
  par(cex.main = 0.65, oma = c(2, 0, 0, 5), mar = c(5, 4.1, 4.1, 5))  # Adjust margins
  coloring <- colorRampPalette(c("blue", "yellow", "red"))
  
  d3heatmap::d3heatmap(as.matrix(mtx), hclust = function(tmp) {
    hclust(tmp, method = input$cmbClustMethod)
  }, colors = coloring(coloring.num), tip_transformation = untransform.method, 
  xaxis_font_size = "10pt", yaxis_font_size = "10pt", xaxis_height = 200, 
  yaxis_height = 200, dendro.rds.path = dend.path)
})
  

output$legendEnrich <- renderPlot({
  if (input$cmbMatrix == "matrix_PVAL.txt") {
    mtx <- get.adjust.matrix()
  } else {
    mtx <- get.matrix()
  }
  coloring <- colorRampPalette(c("blue", "yellow", "red"))
  # Breaks symmetric around 0
  my.breaks <- c(seq(min(mtx), 0, length.out = 10), 0, seq(0, max(mtx), length.out = 10))
  plot(my.breaks, rep(1, length(my.breaks)), col = coloring(length(my.breaks)), 
       pch = 15, cex = 10, main = "Depletion/Enrichment Significance", 
       ylab = "", xlab = "", yaxt = "n", xaxt = "n")
  if (input$cmbMatrix == "matrix_PVAL.txt") {
    axis(side = 1, at = my.breaks, labels = scales::scientific_format(2)(1/10^abs(my.breaks)), las = 1)
  } else {
    axis(side = 1, at = my.breaks, labels = scales::scientific_format(2)(2^my.breaks), las = 1)
  }
})
  

# generate the enrichment table and appends the gf.name information columns
get.enrichment.table <- reactive({
  mtx <- read.csv(paste(get.results.dir(), input$cmbMatrix, sep = ""), sep = "\t")
  selectedFOI <- 1
  selectedFOI <- input$cmbFOI
  if (input$cmbMatrix == "matrix_PVAL.txt") {
    mtx.adjust <- apply(mtx[selectedFOI], 2, 
                        function(x) p.adjust(abs(x), method = input$cmbPvalAdjustMethod)) %>% as.matrix(drop = F)
    pval.sig <- rep("Not significant", nrow(mtx.adjust))
    pval.sig[mtx[, selectedFOI] < 0.05 & mtx[, selectedFOI] > 0] <- "Overrepresented"
    pval.sig[mtx[, selectedFOI] > -0.05 & mtx[, selectedFOI] < 0] <- "Underrepresented"
    
    mtx.table <- cbind(signif(abs(mtx[selectedFOI]), digits = 3), pval.sig, 
                       signif(mtx.adjust, digits = 3))
    colnames(mtx.table) <- c("p.value", "direction", "adj.p.val")
    mtx.table$p.value <- as.numeric(as.character(mtx.table$p.value))
    mtx.table$adj.p.val <- as.numeric(as.character(mtx.table$adj.p.val))
  } else {
    # for odds ratio
    or.sig <- rep("Not significant", nrow(mtx))
    or.sig[mtx[, selectedFOI] > 1] <- "Enriched"
    or.sig[mtx[, selectedFOI] < 1] <- "Depleted"
    mtx.table <- cbind(scales::scientific_format(3)(abs(mtx[selectedFOI])), or.sig)
    colnames(mtx.table) <- c("odds.ratio", "direction")
    mtx.table$odds.ratio <- as.numeric(as.character(mtx.table$odds.ratio))
  }
  # make the GF.name a column instead of just the rowname so we can left_join
  mtx.table = cbind(GF.Name = rownames(mtx.table), mtx.table)  
  rownames(mtx.table) <- NULL
  gfano <- gfAnnot[, (names(gfAnnot) %in% c("file_name", "cell", "cell_desc", 
                                            "factor", "factor_desc", "source", "source_desc"))]
  mtx.table <- dplyr::left_join(mtx.table, gfano, by = c(GF.Name = "file_name"))
  colnames(mtx.table)[1] <- "epigenomic_name"
  return(mtx.table)
})


output$tblEnrichment <- renderDataTable({
  withProgress({
    table.enrich <- get.enrichment.table()
    # shortens descriptions 
    table.enrich$cell_desc <- gr_trimnames(table.enrich$cell_desc, num.char = num.char)
    table.enrich$factor_desc <- gr_trimnames(table.enrich$factor_desc, num.char = num.char)
    table.enrich$source_desc <- gr_trimnames(table.enrich$source_desc, num.char = num.char)
    table.enrich
  }, message = "Loading enrichment table", value = 1)
}, 
options = list(lengthMenu = list(c(10, 50, 100, -1), c("10", "50", "100", "All")), pageLength = 10))


output$downloadEnrichTable <- downloadHandler(filename = function() {
  return("Enrichment_table.txt")
}, content = function(file) {
  write.table(x = get.enrichment.table(), file = file, sep = "\t", quote = F, 
              row.names = F)
})
  
get.annotation.table <- reactive({
  withProgress({
    validate(need(try(mtx <- read.table(paste(get.results.dir(), "annotations/", 
                                              input$cmbAnnotation, sep = ""), header = T)), 
                  "Annotation results not available."))
    mtx
  }, message = "Loading Annotation table", value = 1)
})
  

output$tblAnnotation <- renderDataTable({
  get.annotation.table()
}, 
options = list(lengthMenu = list(c(10, 50, 100, -1), c("10", "50", "100", "All")), pageLength = 10))

outputOptions(output, "downloadEnrichTable", suspendWhenHidden = FALSE)


## enrichment up and down plots for single column --------------------------------------------------------------
get.barplot.matrix <- reactive({
  # populate the enrichment table combobox
  file.names.enrichment <- tools::file_path_sans_ext(list.files(paste(get.results.dir(), "enrichment/", sep = "")))
  mtx <- gr_load_data(paste(get.results.dir(), input$cmbMatrix, sep = ""))
})

axis_text_scaling = 9;
count_barplot = 20;
tbEnrichUp  <-reactive({
  tb.barplot <- tbl_df(data.frame(get.barplot.matrix())) %>% tibble::rownames_to_column(var='gfs')
  if (input$cmbMatrix == "matrix_PVAL.txt") {
    tb.barplot <-  tb.barplot %>% gather(fois,vals,-gfs,convert = T) %>% filter(fois==input$cmbFOI,vals>0) %>% arrange(desc(vals)) %>% slice(1:count_barplot)
  } else {
    tb.barplot <-  tb.barplot %>% gather(fois,vals,-gfs, convert = T) %>% filter(fois==input$cmbFOI,vals>0) %>% arrange(desc(vals)) %>% slice(1:count_barplot)
  }
  # Fill bar char with 0 data rows up till number of less GFs than count_barplot
  le <- length(tb.barplot$vals)
  num <- count_barplot-le
  for (i in seq_len(length.out = num)){
    tb.barplot <-  bind_rows(tb.barplot,data.frame(gfs = toString(i),fois="",vals = 0))
  }
  return(tb.barplot)
})


output$pltEnrichUp_ui <- renderUI({
  tblEnr <- tbEnrichUp()
  maxGFchar <- max(sapply(tblEnr$gfs,function(x) nchar(x)))
  if (length(tblEnr$vals) == 0){
    data.frame(x = 10, y = 10,text = "") %>% ggvis(~x, ~y,text:=~text) %>% layer_text(
      x = prop("x", ~x, scale = "xcenter"),
      y = prop("y", ~y, scale = "ytop"),
      text:=~text, fontSize := 17, fill:="black", baseline:="middle", align:="center"
    ) %>% bind_shiny("pltEnrichDown","pltEnrichDown_ui")
    return("")
  }
  
  tblEnr$gfs <- structure(1:length(tblEnr$vals), .Label = tblEnr$gfs, class = "factor")
  if (input$cmbMatrix == "matrix_PVAL.txt") {
    tblEnr %>% ggvis(x=~gfs,y=~vals, fill := "#e30000") %>% layer_bars(width=.7) %>% 
      add_axis("x",title="", properties = axis_props(labels = list(angle = -45, align = "right", fontSize = 17))) %>% 
      add_axis("y",title="P-value",title_offset = 50) %>%
      set_options(width = 50*length(tblEnr$vals)+100, height = axis_text_scaling * maxGFchar + 200, resizable=T) %>%
      add_tooltip( function(x) {
        if(is.null(x)) return(NULL)
        paste("Genomic Feature: ", x[1],'<br>',"P-value: ", scales::scientific_format(2)(1/10^abs(x[3])),sep=" ")
      }, "hover") %>% 
      bind_shiny("pltEnrichUp","pltEnrichUp_ui")
    
  } else{
    tblEnr %>% ggvis(x=~gfs,y=~vals, fill := "#e30000") %>% layer_bars(width=.7) %>%
     add_axis("x",title="", properties = axis_props(labels = list(angle = -45, align = "right", fontSize = 17))) %>%
      add_axis("y",title="Odds-ratio") %>%
      set_options(width = 50*length(tblEnr$vals)+100, height = axis_text_scaling * maxGFchar + 200, resizable=T) %>%
      add_tooltip( function(x) {
        if(is.null(x)) return(NULL)
        paste("Genomic Feature: ", x[1],'<br>',"Odds-ratio ", scales::scientific_format(2)(2^x[3]),sep=" ")
      }, "hover") %>% 
      bind_shiny("pltEnrichUp","pltEnrichUp_ui")
  }
  return('')
})



tbEnrichDown <-reactive({
  tb.barplot <- tbl_df(data.frame(get.barplot.matrix())) %>% tibble::rownames_to_column(var='gfs')
  if (input$cmbMatrix == "matrix_PVAL.txt") {
    tb.barplot <-  tb.barplot %>% gather(fois,vals,-gfs, convert = T) %>% filter(fois==input$cmbFOI,vals<0) %>% arrange(desc(abs(vals))) %>% slice(1:count_barplot)
  }else{
    # fix for OR
    tb.barplot <-  tb.barplot %>% gather(fois,vals,-gfs, convert = T) %>% filter(fois==input$cmbFOI,vals<0) %>% arrange(desc(abs(vals))) %>% slice(1:count_barplot)
  }
  # Fill bar char with 0 data rows up till number of less GFs than count_barplot
  le <- length(tb.barplot$vals)
  num <- count_barplot-le
  for (i in seq_len(length.out = num)){
    tb.barplot <-  bind_rows(tb.barplot,data.frame(gfs = toString(i),fois="",vals = 0))
  }
  return(tb.barplot)
})



output$pltEnrichDown_ui <- renderUI({
  tblEnr <- tbEnrichDown()
  maxGFchar <- max(sapply(tblEnr$gfs,function(x) nchar(x)))
  if (length(tblEnr$vals) == 0){
    data.frame(x = 10, y = 10,text = "") %>% ggvis(~x, ~y,text:=~text) %>% layer_text(
      x = prop("x", ~x, scale = "xcenter"),
      y = prop("y", ~y, scale = "ytop"),
      text:=~text, fontSize := 17, fill:="black", baseline:="middle", align:="center"
    ) %>% bind_shiny("pltEnrichDown","pltEnrichDown_ui")
    return("")
  }
  tblEnr$gfs <- structure(1:length(tblEnr$vals), .Label = tblEnr$gfs, class = "factor")
  if (input$cmbMatrix == "matrix_PVAL.txt") {
    tblEnr %>% ggvis(x=~gfs,y=~abs(vals), fill := "#1c9600") %>% layer_bars(width=.7) %>%
      add_axis("x",title="", properties = axis_props(labels = list(angle = -45, align = "right", fontSize = 17))) %>%
      add_axis("y",title="P-value",title_offset = 50) %>%
      set_options(width = 50*length(tblEnr$vals)+100, height = axis_text_scaling * maxGFchar + 200, resizable=T) %>%
      add_tooltip( function(x) {
        if(is.null(x)) return(NULL)
        paste("Genomic Feature: ", x[1],'<br>',"P-value: ", scales::scientific_format(2)(1/10^abs(x[3])),sep=" ")
      }, "hover") %>% 
      bind_shiny("pltEnrichDown","pltEnrichDown_ui")
  } else {
  tblEnr %>% ggvis(x=~gfs,y=~abs(vals), fill := "#1c9600") %>% layer_bars(width=.7) %>%
     add_axis("x",title="", properties = axis_props(labels = list(angle = -45, align = "right", fontSize = 17))) %>%
     add_axis("y",title="Odds-ratio") %>%
      set_options(width = 50*length(tblEnr$vals)+100, height =axis_text_scaling * maxGFchar + 200, resizable=T) %>%
      add_tooltip( function(x) {
        if(is.null(x)) return(NULL)
        paste("Genomic Feature: ", x[1],'<br>',"Odds-ratio ",  scales::scientific_format(2)(2^(-x[3])),sep=" ")
     }, "hover") %>% 
      bind_shiny("pltEnrichDown","pltEnrichDown_ui")
  }
  return('') # line required to prevent "ERROR: cannot coerce type 'closure' to vector of type 'character'" from showing
})


  
# episimilarity
# ---------------------------------------------------------------
get.corr.matrix <- reactive({
  mtx <- gr_load_data(paste(get.results.dir(), input$cmbMatrix, sep = ""), p2z = TRUE)
  # If there are columns with SD=0, add jitter to it. Needed for
  # pair-wise column correlation analysis (regulatory similarity
  # analysis). Only valid if there's more than 1 row
  if (nrow(mtx) > 1) {
    ind <- apply(mtx, 2, function(x) sd(x, na.rm = TRUE)) == 0  # Indexes of such columns
    if (sum(ind) > 0) {
      set.seed(1)
      mtx[, ind] <- jitter(mtx[, ind, drop = FALSE], factor = 0.1)
    }
  }
  mtx <- scale(mtx)
  mtx <- Hmisc::rcorr(as.matrix(mtx), type = input$cmbEpisimCorType)[[1]]
  write.table(x = mtx, file = paste(get.results.dir(), "matrix_CORR.txt", sep = ""))
  mtx
})


get.cor.hclust.dendrogram <- reactive({
  cor.mat <- get.corr.matrix()
  as.dendrogram(hclust(as.dist(1 - cor.mat), method = input$cmbClustMethod))
})


output$heatmapEpisim <- renderD3heatmap({
  withProgress({
    mat <- get.matrix()
    validate(need(ncol(mat) > 2, "Need at least 3 SNPs of interest files to perform clustering."))
    validate(need(nrow(mat) > 4, "Need at least 5 genome features to perform clustering."))
    cor.mat <- get.corr.matrix()
    hclustergram <- get.cor.hclust.dendrogram()
  }, message = "Rendering heatmap", value = 1)
  
  coloring <- colorRampPalette(c("blue", "yellow", "red"))
  # TODO: Colv == Rowv?
  d3heatmap::d3heatmap(as.matrix(cor.mat), Rowv = hclustergram, Colv = "Rowv", 
                       colors = coloring(coloring.num), 
                       dendro.rds.path = paste(get.results.dir(), "heatmap.episim.dend.rds", sep = ""), 
                       xaxis_font_size = "10pt", yaxis_font_size = "10pt", xaxis_height = 200, yaxis_height = 200)
  
})


output$legendEpisim <- renderPlot({
  mat <- get.matrix()
  validate(need(ncol(mat) > 2, ""))
  validate(need(nrow(mat) > 4, ""))
  mtx <- get.corr.matrix()
  coloring <- colorRampPalette(c("blue", "yellow", "red"))
  color.range = seq(min(mtx), max(mtx), (max(mtx) - min(mtx))/coloring.num)
  plot(color.range, rep(1, coloring.num + 1), 
       col = coloring(coloring.num + 1), pch = 15, cex = 10, main = "Correlation Coefficient", 
       ylab = "", xlab = "", yaxt = "n")
})


# this function is cut out from the tblEpigenetics renderer. It is a
# long calculation that is only run when # of clusters changes
calculate.clust <- reactive({
  # this line ensure that dendrogram is redrawn when heatmap is
  cor.mat <- get.corr.matrix()  
  # We need the d3heatmap's dendrogram (which is saved to the .rds file)
  mat <- get.matrix()
  validate(need(ncol(mat) > 2, "Need at least 3 SNPs of interest files to perform clustering."))
  validate(need(nrow(mat) > 4, "Need at least 5 genome features to perform clustering."))
  cor.mat <- get.corr.matrix()
  hclustergram <- get.cor.hclust.dendrogram()
  d3heatmap::d3heatmap(as.matrix(cor.mat), Rowv = hclustergram, Colv = "Rowv", 
                       dendro.rds.path = paste(get.results.dir(), "heatmap.episim.dend.rds", sep = ""), 
                       xaxis_font_size = "10pt", yaxis_font_size = "10pt", 
                       xaxis_height = 200, yaxis_height = 200)
  # Dendrogram file is created by d3Heatmap
  dend = readRDS(file = paste(get.results.dir(), "heatmap.episim.dend.rds", sep = ""))  
  # Empirically set desired numter of clusters
  cl_num <- input$sldEpisimNumClust  
  hcut <- dendextend::heights_per_k.dendrogram(dend)[cl_num]  # extract the height to cut based on # of groups
  # get the cluster labels
  mtx.clust <- dend %>% gr_clusters(height = hcut, minmembers = 3)
  mtx = gr_load_data(paste(get.results.dir(), input$cmbMatrix, sep = ""), p2z = FALSE)  # load the original matrix
  is.OR = T
  if (input$cmbMatrix == "matrix_PVAL.txt") {
    is.OR = F
  }
  mtx.deg <- suppressWarnings(gr_degfs(mtx[, mtx.clust$eset.labels], 
                                       mtx.clust, isOR = is.OR))
  updateSelectInput(session, "cmbEpigenetics", "Select which comparison to show", 
                    choices = names(mtx.deg), selected = names(mtx.deg)[1])
  return(mtx.deg)
})


get.epigenetics.table <- reactive({
  withProgress({
    mtx.deg <- calculate.clust()
  }, message = "Calculating regulatory differences", value = 1)
  mtx.deg
})


output$tblEpigenetics <- renderDataTable({
  mtx <- get.matrix()
  validate(need(ncol(mtx) > 2, "Need at least 3 SNPs of interest files to perform clustering."))
  validate(need(nrow(mtx) > 4, "Need at 5 least genome features to perform clustering."))
  # validate(need(try(get.epigenetics.table()),'Either nothing is
  # significant, or there are too few SNP sets per cluster. Re-run the
  # analysis using more SNP sets, or try a different clustering
  # method.'))
  mtx.deg <- get.epigenetics.table()
  selectedCor = input$cmbEpigenetics
  # validate(need(table.epi != 'Results not ready yet.', 'Results not ready yet.'))
  if (selectedCor == "Results not ready yet." | length(mtx.deg) == 0) {
    return(data.frame(NoResult = ""))
  }
  # convert values to numeric form for sorting purposes
  if (input$cmbMatrix == "matrix_PVAL.txt") {
    for (x in list("adj.p.val", 3, 4)) {
      # scientific format returns a factor, use as.character before
      # as.numeric to ensure that the actual number is converted rather than
      # the factor level mtx.deg[[selectedCor]][[x]] <-
      # as.numeric(as.character(scales::scientific_format(3)(as.numeric(mtx.deg[[selectedCor]][[x]]))))
      mtx.deg[[selectedCor]][[x]] <- signif(as.numeric(mtx.deg[[selectedCor]][[x]], digits = 3))
    }
  } else {
    # mtx.deg[[selectedCor]][['adj.p.val']] <-
    # as.numeric(as.character(scales::scientific_format(3)(as.numeric(mtx.deg[[selectedCor]][['adj.p.val']]))))
    mtx.deg[[selectedCor]][["adj.p.val"]] <- signif(as.numeric(mtx.deg[[selectedCor]][["adj.p.val"]],  digits = 3))
  }
  #colnames(mtx.deg[[selectedCor]])[1] <- "epigenomic_name"
  table.epi <- mtx.deg[[selectedCor]][, !(colnames(mtx.deg[[1]]) %in% 
                                            c("full_path", "URL", "full_description", "category", "category_desc"))]
  table.epi$cell_desc <- gr_trimnames(table.epi$cell_desc, num.char = num.char)
  table.epi$factor_desc <- gr_trimnames(table.epi$factor_desc, num.char = num.char)
  table.epi$source_desc <- gr_trimnames(table.epi$source_desc, num.char = num.char)
  table.epi
}, options = list(lengthMenu = list(c(10, 50, 100, -1), c("10", "50", "100", "All")), pageLength = 10))


output$downloadEpigenetics <- downloadHandler(filename = function() {
  return("Regulatory_table.txt")
}, content = function(file) {
  mtx.deg <- get.epigenetics.table()
  selectedCor = input$cmbEpigenetics
  # validate(need(table.epi != 'Results not ready yet.', 'Results not ready yet.'))
  if (mtx.deg == "Results not ready yet." | length(mtx.deg) == 0) {
    return(data.frame(NoResult = "There is nothing signficant to show"))
  }
  # convert values to numeric form for sorting purposes
  if (input$cmbMatrix == "matrix_PVAL.txt") {
    for (x in list("adj.p.val", 3, 4)) {
      mtx.deg[[selectedCor]][[x]] <- scales::scientific_format(3)(as.numeric(mtx.deg[[selectedCor]][[x]]))
    }
  } else {
    mtx.deg[[selectedCor]][["adj.p.val"]] <- scales::scientific_format(3)(as.numeric(mtx.deg[[selectedCor]][["adj.p.val"]]))
  }
  colnames(mtx.deg[[selectedCor]])[1] <- "epigenomic_name"
  table.epi <- mtx.deg[[selectedCor]][, !(colnames(mtx.deg[[1]]) %in% 
                                            c("full_path", "URL", "full_description", "category", "category_desc"))]
  table.epi
  write.table(x = table.epi, file = file, sep = "\t", quote = F, row.names = F)
})


output$pltDend <- renderPlot({
  mtx <- get.matrix()
  validate(need(ncol(mtx) > 2, ""))
  validate(need(nrow(mtx) > 4, ""))
  cor.mat <- get.corr.matrix()  # this line ensure that dendrogram is redrawn when heatmap is
  hclustergram <- get.cor.hclust.dendrogram()  # ensures that dendrogram is redrawn when hclust method is changed
  # Dendrogram file is created by d3Heatmap
  dend = readRDS(file = paste(get.results.dir(), "heatmap.episim.dend.rds", sep = ""))  
  plot(as.dendrogram(dend, hang = -1))  # Plot dendrogram
  cl_num <- input$sldEpisimNumClust  # Empirically set desired numter of clusters
  cols <- rainbow(cl_num)  # Make colos for clusters
  hcut <- dendextend::heights_per_k.dendrogram(dend)[cl_num]  # extract the height to cut based on # of groups
  # get the cluster labels
  mtx.clust <- validate(need(try(dend %>% 
                                   gr_clusters(height = hcut, minmembers = 3)), 
                             "Try using a lower number of clusters"))
  # Define the clusters by rectangles
  validate(need(try(rect.hclust(as.hclust(dend), k = cl_num, border = cols)), 
                "Try a different clustering method."))
})


get.gr_cellspecific <- reactive({
  withProgress({
    mtx <- gr_load_data(paste(get.results.dir(), "matrix_PVAL.txt", sep = ""))
    cellspec_path = paste(get.results.dir(), "cellspecific.rds", sep = "") # this file save the results so we only have to calculate cellspec one time
    if (file.exists(cellspec_path)[1] == F){
      mtx.CTE <-  gr_cellspecific(mtx, cutoff.pval = 0.05)
      saveRDS(mtx.CTE, file=cellspec_path)
    } else {
      mtx.CTE <- readRDS(cellspec_path)
    }
  }, message = "Loading table", value = 1)
  return(mtx.CTE)
})


output$tblCTEnrichment <- renderDataTable({
  inputtmp <- input$cmbFOI
  mtx <- gr_load_data(paste(get.results.dir(), "matrix_PVAL.txt", sep = ""))
  validate(need(nrow(mtx) > 5, "Insufficient data for performing cell type-specific enrichment analysis"))
  # running function
  mtx.CTE <- get.gr_cellspecific()
  if (is.character(mtx.CTE)) {
    return(data.frame(NoResults = "Insufficient data for performing cell type-specific enrichment analysis"))
  }
  if (is.character(mtx.CTE[[input$cmbFOI]])) {
    return(data.frame(NoResults = "Nothing significant"))
  }
  table.CTE <- data.frame(cell = rownames(mtx.CTE[[input$cmbFOI]]), mtx.CTE[[input$cmbFOI]])
  rownames(table.CTE) <- NULL
  # colnames(table.CTE)[2:5] <- c("p.value", "num_of_tests", "av_pval_cell", "av_pval_tot")
  table.CTE$p.value <- as.numeric(table.CTE$p.value)
  table.CTE$av_pval_cell <- as.numeric(table.CTE$av_pval_cell)
  table.CTE$av_pval_tot <- as.numeric(table.CTE$av_pval_tot)
  table.CTE$cell_desc <- gr_trimnames(table.CTE$cell_desc, num.char = num.char)
  table.CTE
}, options = list(lengthMenu = list(c(10, 50, 100, -1), c("10", "50", "100", "All")), pageLength = 10))


output$downloadCTEnrichment <- downloadHandler(filename = function() {
  return("EnrichmentCT_table.txt")
}, content = function(file) {
  selectedCor <- input$cmbFOI
  mtx <- gr_load_data(paste(get.results.dir(), "matrix_PVAL.txt", sep = ""))
  validate(need(nrow(mtx) > 5, "Insufficient data for performing cell type-specific enrichment analysis"))
  # running function
  mtx.CTE <- gr_cellspecific(mtx)
  if (is.character(mtx.CTE)) {
    table.CTE <- data.frame(NoResults = "Insufficient data for performing cell type-specific enrichment analysis")
  } else if (is.character(mtx.CTE[[selectedCor]])) {
    table.CTE <- data.frame(NoResults = "Nothing significant")
  } else {
    table.CTE <- data.frame(cell = rownames(mtx.CTE[[selectedCor]]), mtx.CTE[[selectedCor]])
    rownames(table.CTE) <- NULL
    # colnames(table.CTE)[2:5] <- c("p.value", "num_of_tests", "av_pval_cell", "av_pval_tot")
  }
  write.table(x = table.CTE, file = file, sep = "\t", quote = F, row.names = F)
})
# --download button code--------------------------------------------------

  
output$downloadEnrichHeatmap <- downloadHandler(filename = function() {
  return("EnrichmentHeatmap.pdf")
}, content = function(file) {
  pdf(file = file)
  par(mfrow = c(3, 3))
  if (input$cmbMatrix == "matrix_PVAL.txt") {
    mtx <- get.adjust.matrix()
  } else {
    mtx <- get.matrix()
  }
  n_limit = 30
  # if n > 100, calculate SD for each row.
  if (nrow(mtx) > n_limit) {
    # calculate SD for each row.
    mtx.sd <- apply(mtx, 1, sd)
    mtx.sd <- data.frame(mtx.sd)
    # sort rows based on SD
    mtx.sd.order <- mtx[order(mtx.sd, decreasing = T), ]
    mtx.sd.order <- mtx.sd.order[1:n_limit, ]
    mtx <- mtx.sd.order
  }
  pdf(file = file, width = 0.25*ncol(mtx)+5, height = 0.25*nrow(mtx)+5)
  par(cex.main = 0.65, oma = c(2, 0, 0, 5), mar = c(5, 4.1, 4.1, 5))  # Adjust margins
  gplots::heatmap.2(as.matrix(mtx), hclust = function(tmp) {
    hclust(tmp, method = input$cmbClustMethod)
  }, trace = "none", density.info = "none", 
  col = colorRampPalette(c("blue", "yellow", "red")), main = "Enrichment Heatmap", cexRow = 0.8, cexCol = 1, 
  margins = c(10, 10), srtRow = 0, srtCol = 45)
  dev.off()
}, contentType = "application/pdf")


output$downloadEpisimHeatmap <- downloadHandler(filename = function() {
  return("EpisimilarityHeatmap.pdf")
}, content = function(file) {
  pdf(file = file, width = 10, height = 10)
  mat <- get.matrix()
  validate(need(ncol(mat) > 2, "Need at least 3 SNPs of interest files to perform clustering."))
  validate(need(nrow(mat) > 4, "Need at least 5 genome features to perform clustering."))
  cor.mat <- get.corr.matrix()
  pdf(file = file, width = 0.25*ncol(cor.mat)+5, height = 0.25*nrow(cor.mat)+5)
  hclustergram <- get.cor.hclust.dendrogram()
  par(cex.main = 0.65, oma = c(2, 0, 0, 5), mar = c(5, 4.1, 4.1, 5))  # Adjust margins
  coloring <- colorRampPalette(c("blue", "yellow", "red"))
  gplots::heatmap.2(as.matrix(cor.mat), Colv = hclustergram, Rowv = hclustergram, 
            trace = "none", density.info = "none", 
            col = colorRampPalette(c("blue", "yellow", "red")), main = "Episimilarity Heatmap", 
            cexRow = 0.8, cexCol = 1, margins = c(10, 10), srtRow = 0, srtCol = 45)
  dev.off()
}, contentType = "application/pdf")


output$downloadAnnotation <- downloadHandler(filename = function() {
  return("Enrichment_table.txt")
}, content = function(file) {
  write.table(x = get.annotation.table(), file = file, sep = "\t", quote = F, row.names = F)
})



output$downloadZIP <- downloadHandler(filename = function() {
  return("genome_runner.zip")
}, content = function(file) {
  # ensure correlation matrix is created
  mtx <- get.matrix()
  if (ncol(mtx) > 1 & nrow(mtx) > 4) {
    get.corr.matrix()
  }
  # Append gfAnnot columns to the end of the PVAL and OR matrix
  mtx <- read.csv(paste(get.results.dir(), "matrix_PVAL.txt", sep = ""), sep = "\t")
  mtx <- data.frame(GF = rownames(mtx), mtx)
  mtx <- dplyr::left_join(mtx, gfAnnot, by = c(GF = "file_name"))
  rownames(mtx) <- mtx$GF
  mtx$GF <- NULL
  write.table(mtx, file = paste(get.results.dir(), "matrix_PVAL_annot.txt", 
                                sep = ""), sep = "\t", quote = FALSE, col.names = NA)
  mtx <- read.csv(paste(get.results.dir(), "matrix_OR.txt", sep = ""), sep = "\t")
  mtx <- data.frame(GF = rownames(mtx), mtx)
  mtx <- dplyr::left_join(mtx, gfAnnot, by = c(GF = "file_name"))
  rownames(mtx) <- mtx$GF
  mtx$GF <- NULL
  write.table(mtx, file = paste(get.results.dir(), "matrix_OR_annot.txt", 
                                sep = ""), sep = "\t", quote = FALSE, col.names = NA)
  mtx.clust = tryCatch({
    calculate.clust()
  }, error = function(e) {
    "Try different clustering settings"
  })
  # save each new data frame as an individual .csv file based on its name
  lapply(1:length(mtx.clust), function(i) 
    write.table(mtx.clust[[i]], file = paste0(get.results.dir(), names(mtx.clust[i]), ".txt"), 
                quote = F, row.names = F, sep = "\t"))
  # zip up text files
  files.txt <- as.character(sapply(c("gr_log.txt", "detailed.txt", "matrix_PVAL_annot.txt", 
                                     "matrix_OR_annot.txt", "matrix_CORR.txt"), function(x) {
                                       paste(get.results.dir(), x, sep = "")
                                     }))
  files.txt <- append(files.txt, paste0(get.results.dir(), names(mtx.clust), ".txt"))
  file.names.annotation <- list.files(paste(get.results.dir(), "annotations", sep = ""), full.names = T)
  anot.path <- paste(get.results.dir(), "annotations.zip", sep = "")
  if (length(file.names.annotation) != 0 & !file.exists(anot.path)) {
    zip(zipfile <- anot.path, files = file.names.annotation, flags = "-j")
  }
  enrich.path <- paste(get.results.dir(), "enrichment.zip", sep = "")
  file.names.enrichment <- list.files(paste(get.results.dir(), "enrichment", sep = ""), full.names = T)
  if (length(file.names.enrichment) != 0 & !file.exists(enrich.path)) {
    zip(zipfile <- enrich.path, files = file.names.enrichment, flags = "-j")
  }
  file.all <- c(enrich.path, anot.path, files.txt)
  file.all <- gsub("//", "/", file.all)
  zip(zipfile <- file, files = file.all, flags = "-j")
}, contentType = "application/zip")

  
  # Create a different UI depending if there are multiple GF in the results
  output$mainpage <-renderUI({
    # manually load matrix since controls are not loaded yet
    validate(need(try(mtx <- gr_load_data(paste(get.results.dir(), "matrix_PVAL.txt",sep=""))),"Error loading files. Either no significant results are available, or data files are corrupted. Please, re-run the analysis using larger number of genome annotation datasets."))
    file.names.annotation <- list.files(paste(get.results.dir(),"annotations/",sep=""))
    single.feature = TRUE;
    if (ncol(mtx)>1 & nrow(mtx)>1){single.feature = FALSE}
    if (single.feature == FALSE){
      tabsetPanel(id="tabsMultiple",
                  tabPanel("Enrichment heatmap",
                           br("Heatmap of the enrichment analysis results. Rows - names of regulatory/epigenomic features, shown as \"cell-factor-source\" scheme. Columns - names of SNP sets. Cells contain enrichment p-values/odds ratios. Blue/red gradient highlights depleted/enriched associations, respectively."),
                           br("Mouse over the heatmap to see numerical values. Click-and-drag to zoom in, single click to reset zoom."),
                           downloadButton('downloadEnrichHeatmap',"Download PDF"),
                           d3heatmapOutput("heatmapEnrich", width = "100%", height = "600px"),
                           plotOutput("legendEnrich",width="300px",height="150px")
                  ), 
                  tabPanel("Enrichment barplot",
                           br("Enrichment of the SNP sets in regulatory/epigenomic features, shown as \"cell-factor-source\" names on the X-axis. The height of the bars corresponds to the significance of enriched (top)/depleted (bottom) associations."),
                           br(),
                           h3("Overrepresented",align="center"),
                           ggvisOutput("pltEnrichUp"),
                           uiOutput("pltEnrichUp_ui"),
                           h3("Underrepresented",align="center"),
                           ggvisOutput("pltEnrichDown"),
                           uiOutput("pltEnrichDown_ui")
                  ),
                  tabPanel("Enrichment tables",
                           br("Enrichment analysis results in text format."),
                           br("Table legend: 'epigenomic_name' - names of regulatory/epigenomic features; 'p.value/adj.p.val' - non-adjusted/adjusted for multiple testing p-value of the enrichments; 'direction' - directionality of the enrichments. The following columns describe cell types, factors and regulatory/epigenomic features."),
                           br("Click on a column header (e.g., 'adj.p.val') to sort the table"),
                           br(),
                           downloadButton('downloadEnrichTable', 'Download table'),
                           br(),br(),
                           DT::dataTableOutput("tblEnrichment")),
                  tabPanel("Regulatory similarity heatmap",
                           fluidPage(
                             fluidRow(
                               column(6,
                                      br("Heatmap of regulatory similarity among SNP sets. Cells show correlation coefficients for each pair-wise correlation of the SNP set-specific regulatory profiles."),
                                      br("Mouse over the heatmap to see numerical values. Click-and-drag to zoom in, single click to reset zoom."),
                                      downloadButton('downloadEpisimHeatmap', 'Download PDF'),
                                      d3heatmapOutput("heatmapEpisim", width = "100%", height = "600px"),
                                      plotOutput("legendEpisim",width="300px",height="150px")
                               ),
                               column(6,
                                      br("Dedrogram of regulatory similarity among SNP sets. Clusters define groups of SNP sets with similar regulatory enrichments."),
                                      br("Adjust the number of clusters to identify regulatory/epigenomic differences among them on the \"Differential regulatory analysis\" tab."),
                                      plotOutput("pltDend",width = "100%", height = "500px")
                               )
                             )
                           )),
                  tabPanel("Differential regulatory analysis",
                           br("Differential regulatory analysis identifies regulatory/epigenomic features differentially enriched between clusters of SNP sets (e.g., \"cX_vs_cY\"). Adjust the number of clusters and other clustering metrics on the \"Regulatory similarity heatmap\" tab, if needed."),
                           br(),
                           selectInput("cmbEpigenetics", "Select which comparison to show", choices = list("Results not ready yet.")),
                           p("Table legend: 'epigenomic_name' - names of regulatory/epigenomic features; 'adj.p.val' - adjusted for multiple testing p-value of the enrichment differences between clusters 'cX' and 'cY'; 'cX'/'cY' - average enrichments (p-values or odds ratios) in clusters 'cX' and 'cY', respectively. The following columns describe regulatory/epigenomic features."),
                           downloadButton('downloadEpigenetics', 'Download table'),
                           br(),br(),
                           DT::dataTableOutput("tblEpigenetics")),
                  if (length(file.names.annotation)>0){
                    tabPanel("Annotation Analysis",
                           br("Annotation analysis tables. Here, the overlap between each SNP in a set (rows) and each regulatory/epigenomic feature analyzed (columns) is shown (0=no overlap, 1=overlap)."),
                           br("If more than a 100 regulatory/epigenomic features were selected, the annotation tables are split into multiple tables, each having 100 columns or less."),
                           br(),
                           downloadButton('downloadAnnotation', 'Download table'),
                           br(),br(),
                           DT::dataTableOutput("tblAnnotation"))
                  }else{
                    conditionalPanel('False',tabPanel("Annotation Analysis")
                    )
                  },
                  tabPanel("Cell-type enrichment analysis",
                           br("Cell-type enrichment analysis detects cell type specificity of enrichments of SNP sets. It compares whether a cell type-specific enrichment is significantly different from the overall enrichment of a SNP set. This analysis requires selection of categories with multiple epigenomic/regulatory features per cell type, e.g., \"Histone\" and/or \"chromStates\"."),
                           br("Table legend: \'cell\' - cell type name; \'p.value\' - significance p-value of the differences between average overall enrichment ('av_pval_tot') and average cell type-specific enrichment ('av_pval_cell'); \'num_of_tests\' - how many cell type-specific enrichment tests were used to calculate average cell type-specific enrichment ('av_pval_cell'); 'cell_desc' - description of cell types."),
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
                           br("Enrichment of the SNP sets in regulatory/epigenomic features, shown as \"cell-factor-source\" names on the X-axis. The height of the bars corresponds to the significance of enriched (top)/depleted (bottom) associations."),
                           br(),
                           h3("Overrepresented",align="center"),
                           ggvisOutput("pltEnrichUp"),
                           uiOutput("pltEnrichUp_ui"),
                           h3("Underrepresented",align="center"),
                           ggvisOutput("pltEnrichDown"),
                           uiOutput("pltEnrichDown_ui")
                  ),
                  tabPanel("Enrichment tables",
                           br("Enrichment analysis results in text format."),
                           br("Table legend: 'epigenomic_name' - names of regulatory/epigenomic features; 'p.value/adj.p.val' - non-adjusted/adjusted for multiple testing p-value of the enrichments; 'direction' - directionality of the enrichments. The following columns describe cell types, factors and regulatory/epigenomic features."),
                           br("Click on a column header (e.g., 'adj.p.val') to sort the table"),
                           br(),
                           downloadButton('downloadEnrichTable', 'Download table'),
                           br(),br(),
                           DT::dataTableOutput("tblEnrichment")),
                  if (length(file.names.annotation)>0){
                    tabPanel("Annotation Analysis",
                             br("Annotation analysis tables. Here, the overlap between each SNP in a set (rows) and each regulatory/epigenomic feature analyzed (columns) is shown (0=no overlap, 1=overlap)."),
                             br("If more than a 100 regulatory/epigenomic features were selected, the annotation tables are split into multiple tables, each having 100 columns or less."),
                             br(),
                             downloadButton('downloadAnnotation', 'Download table'),
                             DT::dataTableOutput("tblAnnotation"))
                  }else{
                    conditionalPanel('False',tabPanel("Annotation Analysis")
                    )
                  },
                  tabPanel("Cell-type enrichment analysis",
                           br("Cell-type enrichment analysis detects cell type specificity of enrichments of SNP sets. It compares whether a cell type-specific enrichment is significantly different from the overall enrichment of a SNP set. This analysis requires selection of categories with multiple epigenomic/regulatory features per cell type, e.g., \"Histone\" and/or \"chromStates\"."),
                           br("Table legend: \'cell\' - cell type name; \'p.value\' - significance p-value of the differences between average overall enrichment ('av_pval_tot') and average cell type-specific enrichment ('av_pval_cell'); \'num_of_tests\' - how many cell type-specific enrichment tests were used to calculate average cell type-specific enrichment ('av_pval_cell'); 'cell_desc' - description of cell types."),
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
    validate(need(try(mtx <- gr_load_data(paste(get.results.dir(), "matrix_PVAL.txt",sep=""))),""))
    mtx.col.names <- colnames(data.frame(mtx)) # R converts '-' to '.' in data frames
    file.names.annotation <- list.files(paste(get.results.dir(),"annotations/",sep=""))
    single.feature = TRUE
    if (ncol(mtx)>1 & nrow(mtx)>1){single.feature = FALSE}
    if (single.feature == FALSE){
      sidebarPanel(width = 4,
                   tags$a(href="http://integrativegenomics.org", img(src='GRLogo.png', width="300px")),
                   h3("Data Settings"),
                   conditionalPanel("input.tabsMultiple != 'Cell-type enrichment analysis' && input.tabsMultiple != 'Download'  && input.tabsMultiple != 'Annotation Analysis'",
                     selectInput("cmbMatrix", label = "Results to visualize", 
                                 choices = list("P-values" = "matrix_PVAL.txt", 
                                                "Odds Ratios" = "matrix_OR.txt")),
                     bsTooltip("cmbMatrix", "Select enrichment significance or effect size to visualize/analyze", placement = "top", trigger = "hover")
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
                   bsTooltip("cmbClustMethod", "Select clustering method", placement = "top", trigger = "hover"),
                   conditionalPanel("input.tabsMultiple == 'Regulatory similarity heatmap'",
                                    selectInput('cmbEpisimCorType',label = "Correlation coefficient type",
                                                choices = list("Pearson's" = "pearson",
                                                               "Spearman's" = "spearman"))
                   ),
                  
                   conditionalPanel("input.tabsMultiple == 'Regulatory similarity heatmap'",
                                    hr(),h3("Regulatory similarity"),
                                    sliderInput("sldEpisimNumClust","Number of clusters",min = 2,max=10,value = 4)
                   ),
                   p("Note: Refresh the page if the application stops responding")
      )
    } else { # this is for a single column result file
      sidebarPanel(tags$a(href="http://integrativegenomics.org", img(src='GRLogo.png', width="300px")),
                  h3("Global Settings"), hr(),
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
                   ),
                   p("Note: Refresh the page if the application stops responding")
      )
    }
  })
  # Show UI
  hide("loading_page")
  show("main_content")
  
})

