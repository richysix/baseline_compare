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
  session$userData[['testing']] <- FALSE
  
  # load samples and counts files
  loadedData <- reactive({
    if (session$userData[['debug']]) {
      print('Function: mergedCounts')
    }
    session$userData[['demo']] <- input$demo_data
    if (session$userData[['testing']]) {
      sample_file <- file.path('data', 'Brd2-samples.txt')
      count_file <- file.path('data', 'Brd2-counts.tsv')
    } else if (session$userData[['demo']]) {
      sample_file <- file.path('data', 'Brd2-samples.txt')
      count_file <- file.path('data', 'Brd2-counts.tsv')
    } else{
      sample_file_info <- input$sample_file
      count_file_info <- input$count_file
      if (!is.null(sample_file_info) &
          !is.null(count_file_info)) {
        sample_file <- sample_file_info$datapath
        count_file <- count_file_info$datapath
      } else{
        return(NULL)
      }
    }
    loaded_data <-
      load_data(sample_file, count_file, Mm_baseline, session)
    return(loaded_data)
  })
  
  exptData <- reactive({
    loaded_data <- loadedData()
    if (is.null(loaded_data)) {
      return(NULL)
    } else {
      return(loaded_data[['expt_data']])
    }
  })
  
  combinedData <- reactive({
    loaded_data <- loadedData()
    if (is.null(loaded_data)) {
      return(NULL)
    } else {
      return(loaded_data[['merged_data']])
    }
  })
  
  # check columns
  
  # compare to baseline data
  # check sample stages and gene ids
  
  # combine data, getting appropriate baseline samples by stage
  ddsPlusBaseline <- reactive({
    combined_data <- combinedData()
    if (is.null(combined_data)) {
      return(NULL)
    } else {
      dds <- DESeqDataSet(combined_data, design = ~ sex + condition)
      return(dds)
    }
  })
  
  # run PCA
  pca_info <- reactiveValues()
  output$pca_progress <- renderText({
    dds <- ddsPlusBaseline()
    if (is.null(dds)) {
    } else {
      if (session$userData[['debug']]) {
        print('Function: run_pca')
        print('Variance Stabilizing Transform begin...')
      }
      # Create a Progress object
      progress <- shiny::Progress$new(session)
      # Make sure it closes when we exit this reactive, even if there's an error
      on.exit(progress$close())
      
      progress$set(message = "Calculating PCA...", 
                   detail = 'This will depend on the number of samples', value = 0.3)
      
      dds_vst <- varianceStabilizingTransformation(dds, blind=TRUE)
      
      progress$set(value = 1)
      if (session$userData[['debug']]) {
        cat('Variance Stabilizing Transform done.\n')
      }
      
      pca_info[['pca']] <- prcomp( t( assay(dds_vst) ) )
      pca_info[['propVarPC']] <- pca$sdev^2 / sum( pca$sdev^2 )
      aload <- abs(pca$rotation)
      pca_info[['propVarRegion']] <- sweep(aload, 2, colSums(aload), "/")
      return("PCA completed")
    }
  })
  
  pca_plot_obj <- reactive({
    pca <- pca_info[['pca']]
    if (is.null(pca)) {
      return(NULL)
    } else {
      plot_data <- data.frame(
        pc1 = pca[['x']][,1],
        pc2 = pca$x[,2],
        shape = factor(c(rep('het', 6), rep('hom', 5), rep('wt', 3), rep('baseline', 111)),
                       levels = c('baseline', 'wt', 'het', 'hom')),
        colour = colData(dds_vst)[['stage']]
      )
      pca_plot <- ggplot(data = plot_data) + 
        geom_point( aes(x = pc1, y = pc2, shape = shape, fill = colour), size = 3) + 
        scale_shape_manual(values = c(21:24)) + 
        guides(fill = guide_legend(override.aes = list(shape = 21)))
      
      return(pca_plot)
    }
  })
  
  # render heatmap plot
  output$pca_plot <- renderPlot({
    return(pca_plot_obj())
  })
  
  # run DESeq2
  
  # overlap the 3 DE lists
  
  # show results in results tab
}
