library(shiny)
library(SummarizedExperiment)
library(DESeq2)
library(ggplot2)
library(shinyMisc)
library(shinyBS)

# source functions
source(file.path('R', 'load_data.R'))
source(file.path('R', 'deseq_functions.R'))

# load baseline data
# loads object called Mm_baseline
load(file.path('data', 'Mm_baseline_data.rda'))

# Server logic
server <- function(input, output, session) {
  # set testing and debugging options
  session$userData[['debug']] <- TRUE
  session$userData[['testing']] <- TRUE
  
  if (session$userData[['testing']]) {
    load('data/test-baseline.rda')
    Mm_baseline <- Mm_baseline_test
  }
  
  # load samples and counts files
  exptData <- reactive({
    if (session$userData[['debug']]) {
      print('Function: exptData')
    }
    session$userData[['demo']] <- input$demo_data
    if (session$userData[['testing']]) {
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
      load_data(sample_file, count_file, session)
    return(loaded_data)
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
  
  # create DESeq2 data sets
  deseqDatasets <- reactive({
    button_value <- input$analyse_data
    if (button_value > 0) {
      expt_data <- isolate(exptData())
      if (is.null(expt_data)) {
        return(NULL)
      } else {
        use_gender <- isolate(input$use_gender)
        if ( use_gender ) {
          gender_column = isolate(input$sex_var)
          groups <- c('sex')
        } else {
          gender_column = NULL
          groups <- NULL
        }
        # experiment data only
        dds <- 
          create_new_DESeq2DataSet(expt_data, baseline_data = NULL, gender_column = gender_column,
                                   groups = groups, condition_column = isolate(input$condition_var), 
                                   session_obj = session )
        
        # experiment data plus stage matched baseline samples
        expt_plus_baseline_dds <-
          withCallingHandlers(
            create_new_DESeq2DataSet(expt_data, baseline_data = Mm_baseline, gender_column = gender_column,
                                     groups = groups, condition_column = isolate(input$condition_var), 
                                     match_stages = TRUE, session_obj = session ),
            error = function(e){ stop(e) },
            warning = function(w){
              print(conditionMessage(w))
              createAlert(session, anchorId = 'input_file_alert', 
                          content = conditionMessage(w), style = 'warning')
              invokeRestart("muffleWarning")
            })
        
        # experiment data plus all baseline samples
        expt_plus_all_baseline_dds <-
          suppressWarnings(
            create_new_DESeq2DataSet(expt_data, baseline_data = Mm_baseline, gender_column = gender_column,
                                     groups = groups, condition_column = isolate(input$condition_var), 
                                     match_stages = FALSE, session_obj = session ))

        # experiment data plus stage matched baseline samples
        # design formula includes stage
        groups <- c(groups, 'stage')
        expt_plus_baseline_with_stage_dds <-
          suppressWarnings(
            create_new_DESeq2DataSet(expt_data, baseline_data = Mm_baseline, gender_column = gender_column,
                                     groups = groups, condition_column = isolate(input$condition_var), 
                                     match_stages = TRUE, session_obj = session ))

        return(
          list(
            dds = dds,
            expt_plus_baseline_dds = expt_plus_baseline_dds,
            expt_plus_all_baseline_dds = expt_plus_all_baseline_dds,
            expt_plus_baseline_with_stage_dds = expt_plus_baseline_with_stage_dds
          )
        )
      }
    }
  })
  
  # run PCA with just match baseline samples for speed
  pca_info <- reactiveValues()
  output$pca_progress <- renderText({
    expt_data <- exptData()
    deseq_datasets <- deseqDatasets()
    if (is.null(expt_data)) {
      return('Loading Data ...')
    } else {
      if (is.null(deseq_datasets)) {
        return('Experiment Data loaded.')
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
        
        dds_vst <- varianceStabilizingTransformation(deseq_datasets[['expt_plus_baseline_dds']], blind=TRUE)
        
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
        return("PCA completed.")
      }
    }
  })

  pca_plot_obj <- reactive({
    pca <- pca_info[['pca']]
    if (is.null(pca)) {
      return(NULL)
    } else {
      dds_vst <- pca_info[['dds_vst']]
      plot_data <- data.frame(
        PC1 = pca[['x']][,1],
        PC2 = pca[['x']][,2],
        shape = colData(dds_vst)[['condition']],
        colour = colData(dds_vst)[['stage']]
      )
      pca_plot <- ggplot(data = plot_data) +
        geom_point( aes(x = PC1, y = PC2, shape = shape, fill = colour), size = 3) +
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
