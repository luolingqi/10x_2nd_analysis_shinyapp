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
# UI starts here!
###################################################################################################
ui <- fluidPage(theme = shinytheme("cosmo"),
                
  #titlePanel(verbatimTextOutput("The shiny app for 10X genomics secondary analysis \n\n ")),
  #shinyjs::inlineCSS(list(body = "color: red;font-size: 50px;font-style: italic;font-weight: bold;")),
  titlePanel(verbatimTextOutput("title")),
  tags$head(tags$style("#title{color: red;
                               font-size: 50px;
                               #font-style: italic;
                               font-weight: bold;
                              }"
                      )
            ),
  tabsetPanel(
    tabPanel("Quality Control",fluid = TRUE,
      sidebarLayout(
        sidebarPanel(h2("Section A: Quality Control Analysis"),
                     width = 12
                    ),
        mainPanel(
          
          do.call(tabsetPanel, c(id = "q", 
                                 lapply(2:3, function(i) 
                                   {
                                     tabPanel(tabPanel_name[i],
                                              fluidRow(
                                                column(8,
                                                  h3(tabPanel_title[i]),
                                                  plotOutput(tabPanel_name[i]),
                                                  p(tabPanel_legend[i]),
                                                  downloadButton(paste0("bn_",tabPanel_name[i]),"Download")
                                                ),
                                                column(4,
                                                    textOutput(paste0("text_",tabPanel_name[i]))
                                                )
                                             )
                                     )
                                    }
                                   )
                                )
                  )

        )
      ) # closing sidebarLayout
    ), # closing QC tabPanel
    tabPanel("Exploratory Analysis", fluid = TRUE,
             sidebarLayout(
               sidebarPanel(
                 h2("Section B: Exploratory Analysis"),
                 width = 12
               ),
               mainPanel(
                 do.call(tabsetPanel, c(id = "e", lapply(4:13, function(i) {
                   tabPanel(tabPanel_name[i],
                            h3(tabPanel_title[i]),
                            plotOutput(tabPanel_name[i]),
                            p(tabPanel_legend[i]),
                            downloadButton(paste0("bn_",tabPanel_name[i]),"Download")
                   )
                   
                 })
                 )
                 )
                 
                 
                 
               ) # end of mainPanel
             ) # end of sidebarLayout
    )
  ) # closing the outer tabsetPanel
)










