library(shiny)
library(shinythemes)
library(shinyjs)
#library(ggplot2)

library("DropletUtils")
library(scater)
library(scran)
library(EnsDb.Hsapiens.v86)


###################################################################################################
# Prepare text content
###################################################################################################

tabPanel_name <- c("summary_MT", 
                   "umi_cell_ranking",
                   "umi_feature_count_distribution",
                   "average_gene_expression",
                   "most_highly_expressed_genes",
                   "total_count_vs_size_factor",
                   "mean_variance_plot",
                   "expression_distribution_for_top_10_largest_biological_components",
                   "Variance_explained_by_PCs",
                   "Pairwise_plot_for_top3_PCs",
                   "tSNE_plot_from_denoised_PCs",
                   "heatmap_of_clusters",
                   "examine_marker_genes")
tabPanel_title <- c("Summary of Mitochondria genes", 
                    "Total UMI count vs cell ranking",
                    "Plots for cell quality",
                    "Average gene expression",
                    "Most highly expressed genes",
                    "Total UMI count vs size factor",
                    "Mean Variance Plot",
                    "Expression distribution for top 10 largest biological components",
                    "Variance explained by PCs",
                    "Pairwise plot for the top3 PCs",
                    "tSNE plot from denoised PCs",
                    "Heatmap of expression by clusters",
                    "Examine the marker genes")
tabPanel_legend <- c("Table 1. This is a brief summary table, showing the number of Mitochondria genes",
                     "Figure 1. This plot ranks the cells according to the total number of UMI counts. The knee point roughly estimates the cell number",
                     "Figure 2. This plot demonstrates the overall cell quality of the sample by displaying the distribution of
                      1) total intracellular UMI count and   2) total features (genes) count",
                     "Figure 3. This plot demonstrate the average gene expression level as a QC",
                     "Figure 4. This plot shows mostly expressed genes across all the cells",
                     "Figure 5. Plot of total UMI count vs size factor",
                     "Figure 6. Mean variance plot across all the genes",
                     "Figure 7. Top 10 genes with biggest variation",
                     "Figure 8. Variance explained by the principle components",
                     "Figure 9. Pairwise plots for the top3 principle components",
                     "Figure 10. tSNE plot from denoised PCs",
                     "Figure 11. Heatmap of expression by clusters",
                     "Figure 12. Heatmap of expression for marker genes")





###################################################################################################
# Server starts here!
###################################################################################################
server <- function(input, output) {
  
  ######## Download Handler ###############################################
  # Save the file to local server folder and file copy to customer computer
  ##########################################################################
  Downloadfile <- function(graph,plotFunc) {
    # Save the file locally on the Shiny server folder
    pdf(paste0(graph,".pdf"))
    #eval(paste0("plot_",graph,"()"))
    plotFunc()
    dev.off()
        
    output[[paste0("bn_",graph)]] <- downloadHandler(
      filename = function() {
        paste0(graph,".pdf")
      },
      content = function(file) {
        file.copy(paste0(graph,".pdf"), file, overwrite = TRUE)
      }
    )
  }
  
  output$title <- renderText({
    paste("Secondary analysis for 10X genomics scRNASeq data ","\n", sep = "\n")
  })

  # output$hist <- renderPlot({
  #   hist(rnorm(input$num))
  # })
  matrix_path <- "./GRCh38"
  sce <- read10xCounts(matrix_path, col.names = TRUE)
  #

  # Reference: https://academic.oup.com/bioinformatics/article/33/8/1179/2907823
  #
  # uniquify the feature name to remove duplicate and missing values
  rownames(sce) <- uniquifyFeatureNames(rowData(sce)$ID, rowData(sce)$Symbol)
  # ADD chromosome to rowData of sce (Human Genome only!)

  location <- mapIds(EnsDb.Hsapiens.v86, keys = rowData(sce)$ID, column = "SEQNAME", keytype = "GENEID")
  rowData(sce)$CHR <- location
  summary_MT<-as.data.frame(table(location=="MT",useNA = "ifany"))
  colnames(summary_MT) <- c("MT_or_not","Freq")
  
  #Dynamically insert a tabPanel to an existing tabsetPanel
  ##########################################################################
  insertTab(inputId = "q",
            tabPanel("Summary (Mitochondria Genes)", 
                             h3("Summary of Mitochondria genes"),
                             tableOutput("summary_MT"),
                             p("Table 1. This is a brief summary table, showing the number of Mitochondria genes")
                             ),
            target = "umi_cell_ranking",
            position = "before")
  
  output$summary_MT <- renderTable(summary_MT)
  ##########################################################################
  ## 2 Call cells from empty droplets
  bcrank <- barcodeRanks(counts(sce))
  # Only showing unique points for plotting speed.
  uniq <- !duplicated(bcrank$rank)
  #
  plot_umi_cell_ranking <- function(){
    ## Plot total UMI count for each barcode in the dataset against its rank (in decreasing order of total counts)
    plot(bcrank$rank[uniq], bcrank$total[uniq], log="xy",
         xlab="Rank", ylab="Total UMI count", cex.lab=1.2, cex = 0.5, xaxt='n',yaxt='n')
    axis(1, at = c(1,10,100,1000,2000,10000,50000,100000), las=3)
    axis(2, at = c(1,10,100,1000,10000))
    #
    abline(h=bcrank$inflection, col="darkgreen", lty=2)
    abline(h=bcrank$knee, col="dodgerblue", lty=2)
    #
    legend("bottomleft", legend=c("Inflection", "Knee"),
           col=c("darkgreen", "dodgerblue"), lty=2, cex=1.2)
  }
  output$umi_cell_ranking <- renderPlot({
    plot_umi_cell_ranking()
    })
  # try the textOut
  #output$text_umi_cell_ranking <- renderText({"hello how are you?"})
  
  Downloadfile("umi_cell_ranking", plot_umi_cell_ranking)
  
  ## 3 Predict/extract the non-empty drops out of ambient background
  set.seed(100)
  e.out <- emptyDrops(counts(sce))
  sum(e.out$FDR <= 0.01, na.rm=TRUE) # 2004 cells have the expression profile significantly different from the ambient pool
  print("Summary table for empty drop analysis at the FDR cutoff 0.01")
  #table(Sig=e.out$FDR <= 0.01, Limited=e.out$Limited)
  summary_empty_drop<-as.data.frame(table(Sig=e.out$FDR <= 0.01, Limited=e.out$Limited))
  #colnames(summary_MT) <- c("MT_or_not","Freq")
  #Dynamically insert a tabPanel to an existing tabsetPanel
  ##########################################################################
  insertTab(inputId = "q",
            tabPanel("Prediction of empty drops", 
                     h3("Prediction of empty drops"),
                     tableOutput("summary_empty_drop"),
                     p("Table 2. This is a brief summary table of empty-drop prediction")
            ),
            target = "umi_feature_count_distribution",
            position = "before")
  
  output$summary_empty_drop <- renderTable(summary_empty_drop)
  ##########################################################################
  # using which() to automatically remove NAs.
  sce <- sce[,which(e.out$FDR <= 0.01)]
  #
  ## 4 statement about the quality of the cells
  ## in order to identify and remove the droplets containing damaged or dying cells before downstream analysis
  #sce <- calculateQCMetrics(sce)
  sce <- calculateQCMetrics(sce, feature_controls=list(Mito=which(location=="MT")))

  plot_umi_feature_count_distribution <- function() {
    par(mfrow=c(3,1))
    hist(sce$log10_total_counts, breaks=20, col="grey80",
         xlab="Log-total UMI count")
    hist(sce$log10_total_features_by_counts, breaks=20, col="grey80",
         xlab="Log-total number of expressed features")

  }
  output$umi_feature_count_distribution <- renderPlot({
    plot_umi_feature_count_distribution()
    })
  Downloadfile("umi_feature_count_distribution", plot_umi_feature_count_distribution)
  #
  # remove cells with large mitochondrial proportions, representing "CELL DAMAGE!!!"
  high.mito <- isOutlier(sce$pct_counts_Mito, nmads=3, type="higher")
  sce <- sce[,!high.mito]
  summary(high.mito)
  summary_high_mito <- as.data.frame(table(high.mito))
  
  #Dynamically insert a tabPanel to an existing tabsetPanel
  #################################################################
  insertTab(inputId = "q",
            tabPanel("Cells with high mitochondrial portion", 
                     h3("Cells with high mitochondrial portion"),
                     tableOutput("summary_high_mito"),
                     p("Table 3. This is a brief summary table, showing cells with high mito portion")
            ),
            target = "umi_feature_count_distribution",
            position = "after")
  
  output$summary_high_mito <- renderTable(summary_high_mito)
  #################################################################
  ## 5 Examine gene expression after normalization of library size/ size factor
  ave <- calcAverage(sce)
  rowData(sce)$AveCount <- ave
  
  ## Histogram of the log10-average counts for each gene in the dataset
  plot_average_gene_expression <- function() {
    hist(log10(ave), col = "grey80")
  }
  output$average_gene_expression <- renderPlot({
    plot_average_gene_expression()
  })
  # Set up the download button handler by directly copying the file on the server
  Downloadfile("average_gene_expression", plot_average_gene_expression)
  
  ## Plot of the 50 most highly expressed genes
  plot_most_highly_expressed_genes <- function() {
    plotHighestExprs(sce)
  }
  output$most_highly_expressed_genes <- renderPlot({
    plot_most_highly_expressed_genes()
  })
  # Set up the download button handler by directly copying the file on the server
  Downloadfile("most_highly_expressed_genes", plot_most_highly_expressed_genes)
  
  # output$QC_plots <- renderPlot({
  #   plotQC(sce)
  # })
  
  ## 6 Normalize for cell-specific biases
  #  quick preclustering
  
  clusters <- quickCluster(sce, method="igraph", min.mean=0.1, irlba.args=list(maxit=1000)) # for convergence.
  print("A quick pre-clustering to avoid pooling cells that are very different")
  table(clusters)
  summary_quick_clustering <- as.data.frame(table(clusters))
  #Dynamically insert a tabPanel to an existing tabsetPanel
  #################################################################
  insertTab(inputId = "e",
            tabPanel("A quick pre-clustering", 
                     h3("A quick pre-clustering"),
                     tableOutput("summary_quick_clustering"),
                     p("Table 4. This is a quick pre-clustering to avoid pooling cells that are very different")
            ),
            target = "total_count_vs_size_factor",
            position = "before")
  
  output$summary_quick_clustering <- renderTable(summary_quick_clustering)
  #################################################################
  #
  #  Deconvolution method to compute size factors for all cells
  #  Publication Ref: https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0947-7
  sce <- computeSumFactors(sce, min.mean=0.1, cluster=clusters)
  #  Summary of size factor
  # summary(sizeFactors(sce))
  # summary_size_factor <- as.data.frame(summary(sizeFactors(sce)))
  # #Dynamically insert a tabPanel to an existing tabsetPanel
  # #################################################################
  # insertTab(inputId = "e",
  #           tabPanel("Size factor distribution", 
  #                    h3("Size factor distribution"),
  #                    tableOutput("summary_size_factor"),
  #                    p("Table 5. This table gives the idea of size factor distribution")
  #           ),
  #           target = "mean_variance_plot",
  #           position = "before")
  # 
  # output$summary_size_factor <- renderTable(summary_size_factor)
  # #################################################################
  
  #
  # Plot the total count against size factor
  plot_total_count_vs_size_factor <- function() {
    plot(sce$total_counts, sizeFactors(sce), log="xy")
  }
  output$total_count_vs_size_factor <- renderPlot({
    plot_total_count_vs_size_factor()
  })
  Downloadfile("total_count_vs_size_factor", plot_total_count_vs_size_factor)
  
  sce <- normalize(sce)
  #
  ## 7 Modeling the mean-variance trend
  #xxxxxxxx  There is sth wrong with the function 'makeTechTrend',
  #xxxxxxxx     Error in seq.default(from = 0, to = upper.value, length.out = 100) :
  #xxxxxxxx     'to' must be a finite number
  #
  new.trend <- makeTechTrend(x=sce)
  #
  fit <- trendVar(sce, use.spikes=FALSE, loess.args=list(span=0.05))
  plot_mean_variance_plot <- function() {
    plot(fit$mean, fit$var, pch=16)
    curve(fit$trend(x), col="dodgerblue", add=TRUE)
    curve(new.trend(x), col="red", add=TRUE)
  }
  output$mean_variance_plot <- renderPlot({
    plot_mean_variance_plot()
  })
  Downloadfile("mean_variance_plot", plot_mean_variance_plot)
  
  #
  fit0 <- fit
  fit$trend <- new.trend
  dec <- decomposeVar(fit=fit)
  top.dec <- dec[order(dec$bio, decreasing=TRUE),]
  head(top.dec)
  #
  plot_expression_distribution_for_top_10_largest_biological_components <- function() {
    plotExpression(sce, features=rownames(top.dec)[1:10])
  }
  output$expression_distribution_for_top_10_largest_biological_components <- renderPlot({
    plot_expression_distribution_for_top_10_largest_biological_components()
  })
  Downloadfile("expression_distribution_for_top_10_largest_biological_components", plot_expression_distribution_for_top_10_largest_biological_components)
  
  # 8 Dimensionality reduction
  sce <- denoisePCA(sce, technical=new.trend, approx=TRUE)
  print("The number of dimensions to retain after PCA")
  ncol(reducedDim(sce, "PCA"))
  
  plot_Variance_explained_by_PCs <- function() {
    plot(attr(reducedDim(sce), "percentVar"), xlab="PC",
         ylab="Proportion of variance explained")
    abline(v=ncol(reducedDim(sce, "PCA")), lty=2, col="red")
  }
  output$Variance_explained_by_PCs <- renderPlot({
    plot_Variance_explained_by_PCs()
  })
  Downloadfile("Variance_explained_by_PCs", plot_Variance_explained_by_PCs)
  
  plot_Pairwise_plot_for_top3_PCs <- function() {
    plotPCA(sce, ncomponents=3, colour_by="log10_total_features_by_counts")
  }
  output$Pairwise_plot_for_top3_PCs <- renderPlot({
    plot_Pairwise_plot_for_top3_PCs()
  })
  Downloadfile("Pairwise_plot_for_top3_PCs", plot_Pairwise_plot_for_top3_PCs)
  
  sce <- runTSNE(sce, use_dimred="PCA", perplexity=30, rand_seed=100)
  
  plot_tSNE_plot_from_denoised_PCs <- function() {
    plotTSNE(sce, colour_by="log10_total_features_by_counts")
  }
  output$tSNE_plot_from_denoised_PCs <- renderPlot({
    plot_tSNE_plot_from_denoised_PCs()
  })
  Downloadfile("tSNE_plot_from_denoised_PCs", plot_tSNE_plot_from_denoised_PCs)
  
  # 9 Clustering with graph-based methods
  # We build a shared nearest neighbour graph (Xu and Su 2015) and use the Walktrap algorithm to identify clusters.
  snn.gr <- buildSNNGraph(sce, use.dimred="PCA")
  clusters <- igraph::cluster_walktrap(snn.gr)
  sce$Cluster <- factor(clusters$membership)
  print("summary table of a graph-based clustering")
  table(sce$Cluster)
  
  summary_graph_cluster <- as.data.frame(table(sce$Cluster))
  #Dynamically insert a tabPanel to an existing tabsetPanel
  #################################################################
  insertTab(inputId = "e",
            tabPanel("Graph based clustering", 
                     h3("Graph based clustering"),
                     tableOutput("summary_graph_cluster"),
                     p("Table 6. This is the summary result of graph base clustering")
            ),
            target = "heatmap_of_clusters",
            position = "before")
  
  output$summary_graph_cluster <- renderTable(summary_graph_cluster)
  #################################################################
  cluster.mod <- clusterModularity(snn.gr, sce$Cluster, get.values=TRUE)
  log.ratio <- log2(cluster.mod$observed/cluster.mod$expected + 1)
  
  library(pheatmap)
  plot_heatmap_of_clusters <- function() {
    pheatmap(log.ratio, cluster_rows=FALSE, cluster_cols=FALSE, color=colorRampPalette(c("white", "blue"))(100))
  }
  output$heatmap_of_clusters <- renderPlot({
    plot_heatmap_of_clusters()
  })
  Downloadfile("heatmap_of_clusters", plot_heatmap_of_clusters)
  
  # output$tSNE_plot_from_denoised_PCs_1 <- renderPlot({
  #   plotTSNE(sce, colour_by="Cluster")
  # })
  
  ## 10 Marker gene detection
  markers <- findMarkers(sce, clusters=sce$Cluster, direction="up")
  
  # examine the markers for cluster 1
  print("Examine the markers for cluster 1")
  marker.set <- markers[["1"]]
  head(marker.set[,1:8], 10) # only first 8 columns, for brevity
  
  chosen <- rownames(marker.set)[marker.set$Top <= 10]
  
  plot_examine_marker_genes <- function() {
    plotHeatmap(sce, features=chosen, exprs_values="logcounts",
                zlim=5, center=TRUE, symmetric=TRUE, cluster_cols=FALSE,
                colour_columns_by="Cluster", columns=order(sce$Cluster))
  }
  output$examine_marker_genes <- renderPlot({
    plot_examine_marker_genes()
  })
  Downloadfile("examine_marker_genes", plot_examine_marker_genes)

}
