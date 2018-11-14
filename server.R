library(shiny)
library(SummarizedExperiment)
library(DESeq2)
library(ggplot2)

# source functions
source(file.path('R', 'load_data.R'))

# load baseline data
# loads object called baseline_data
load(file.path('data', 'Mm_baseline_data.rda'))

# Server logic
server <- function(input, output, session) {
  # set testing and debugging options
  session$userData[['debug']] <- TRUE
  session$userData[['testing']] <- TRUE
  
  # load samples and counts files
  combinedData <- reactive({
    if (session$userData[['debug']]) {
      print('Function: mergedCounts')
    }
    if (session$userData[['testing']]) {
      sample_file <- file.path('data', 'Brd2-samples.txt')
      count_file <- file.path('data', 'Brd2-counts.tsv')
      combined_data <-
        load_data(sample_file, count_file, Mm_baseline, session)
    } else{
      # sample_file_info <- input$sample_file
      # count_file_info <- input$count_file
      # if (!is.null(sampleFileInfo) &
      #     !is.null(countFileInfo)) {
      #   exptData <-
      #     load_data(sampleFileInfo$datapath,
      #               countFileInfo$datapath,
      #               input$dataType,
      #               session)
      # } else if (!is.null(dataFileInfo)) {
      #   load(dataFileInfo$datapath)
      # }
      # else{
      #   return(NULL)
      # }
    }
    return(combined_data)
  })
  
  # check columns
  
  # compare to baseline data
  # check sample stages and gene ids
  
  # combine data, getting appropriate baseline samples by stage
  dds_plus_baseline <- reactive({
    combined_data <- combinedData()
    dds <- DESeqDataSet(combined_data, design = ~ sex + condition)
    return(dds)
  })
  
  # run PCA
  pca_plot_obj <- reactive({
    dds <- dds_plus_baseline()
    if (session$userData[['debug']]) {
      print('Function: pca_plot_obj')
      print('Variance Stabilizing Transform begin...')
    }
    
    dds_vst <- varianceStabilizingTransformation(dds, blind=TRUE)
    
    if (session$userData[['debug']]) {
      cat('Variance Stabilizing Transform done.\n')
    }
    
    pca <- prcomp( t( assay(dds_vst) ) )
    propVarPC <- pca$sdev^2 / sum( pca$sdev^2 )
    aload <- abs(pca$rotation)
    propVarRegion <- sweep(aload, 2, colSums(aload), "/")
    plot_data <- data.frame(
      pc1 = pca$x[,1],
      pc2 = pca$x[,2],
      shape = factor(c(rep('het', 6), rep('hom', 5), rep('wt', 3), rep('baseline', 111)),
                     levels = c('baseline', 'wt', 'het', 'hom')),
      colour = colData(dds_vst)[['stage']]
    )
    pca_plot <- ggplot(data = plot_data) + 
      geom_point( aes(x = pc1, y = pc2, shape = shape, fill = colour)) + 
      scale_shape_manual(values = c(21:24)) + 
      guides(fill = guide_legend(override.aes = list(shape = 21)))
    
    return(pca_plot)
  })
  
  # render heatmap plot
  output$pca_plot <- renderPlot({
    return(pca_plot_obj())
  })
  
  # run DESeq2
  
  # overlap the 3 DE lists
  
  # show results in results tab
}
