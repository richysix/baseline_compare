library(shiny)
library(SummarizedExperiment)
library(DESeq2)
library(ggplot2)
library(shinyMisc)

# source functions
source(file.path('R', 'load_data.R'))

# load baseline data
# loads object called Mm_baseline
load(file.path('data', 'Mm_baseline_data.rda'))

# Server logic
server <- function(input, output, session) {
  # set testing and debugging options
  session$userData[['debug']] <- TRUE
  session$userData[['testing']] <- TRUE
  
  # load samples and counts files
  loadedData <- reactive({
    if (session$userData[['debug']]) {
      print('Function: loadedData')
    }
    session$userData[['demo']] <- input$demo_data
    if (session$userData[['testing']]) {
      load('data/test-baseline.rda')
      Mm_baseline <- Mm_baseline_test
      sample_file <- file.path('data', 'test-brd2-samples.txt')
      count_file <- file.path('data', 'test-brd2-counts.tsv')
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
  
  factors_in_data <- reactive({
    expt_data <- exptData()
    if (is.null(expt_data)) {
      return(NULL)
    } else {
      return( shinyMisc::factors_in_data(colData(expt_data)) )
    }
  })
  
  # Change UI options based on input data
  # Sex column
  observe({
    if (session$userData[['debug']]) {
      cat("Function: UI observer - Col Names\n")
    }
    if (!input$use_gender) {
      updateRadioButtons(session, "sex_var",
                         choices = list("None"),
                         selected = "None"
      )
    }
    if (!is.null(factors_in_data()) & input$use_gender) {
        # sex options
        sex_options <- as.list(factors_in_data())
        names(sex_options) <- factors_in_data()
        if ( any(sex_options == 'sex') ) {
          selected_sex_option <- 'sex'
        } else {
          selected_sex_option <- sex_options[[1]]
        }
        updateRadioButtons(session, "sex_var",
                           choices = sex_options,
                           selected = selected_sex_option
        )
    }
    if (!is.null(factors_in_data())) {
        # condition options
        condition_options <- as.list(factors_in_data())
        names(condition_options) <- factors_in_data()
        updateRadioButtons(session, "condition_var",
                           choices = condition_options,
                           selected = condition_options[[1]]
        )
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

      pca <- prcomp( t( assay(dds_vst) ) )
      pca_info[['pca']] <- pca
      pca_info[['propVarPC']] <- pca$sdev^2 / sum( pca$sdev^2 )
      aload <- abs(pca$rotation)
      pca_info[['propVarRegion']] <- sweep(aload, 2, colSums(aload), "/")
      pca_info[['dds_vst']] <- dds_vst
      return("PCA completed")
    }
  })

  pca_plot_obj <- reactive({
    pca <- pca_info[['pca']]
    dds_vst <- pca_info[['dds_vst']]
    if (is.null(pca)) {
      return(NULL)
    } else {
      plot_data <- data.frame(
        pc1 = pca[['x']][,1],
        pc2 = pca$x[,2],
        shape = colData(dds_vst)[['condition']],
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
