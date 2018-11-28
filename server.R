library(shiny)
library(shinyBS)
library(DT)
library(SummarizedExperiment)
library(DESeq2)
library(ggplot2)
library(scales)
library(shinyMisc)
library(biovisr)

# source functions
source(file.path('R', 'load_data.R'))
source(file.path('R', 'deseq_functions.R'))
source(file.path('R', 'helper_functions.R'))

# Server logic
server <- function(input, output, session) {
  # load baseline data
  # loads object called Mm_baseline
  load(file.path('data', 'Mm_baseline_data.rda'))
  
  # set testing and debugging options
  session$userData[['debug']] <- TRUE
  session$userData[['testing']] <- TRUE
  
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
        condition_column <- isolate(input$condition_var)
        if ( use_gender ) {
          gender_column <- isolate(input$sex_var)
          groups <- c('sex')
        } else {
          gender_column = NULL
          groups <- NULL
        }
        # experiment data only
        expt_only_dds <- 
          create_new_DESeq2DataSet(expt_data, baseline_data = NULL, gender_column = gender_column,
                                   groups = groups, condition_column = condition_column, 
                                   session_obj = session )
        
        # experiment data plus stage matched baseline samples
        expt_plus_baseline_dds <-
          withCallingHandlers(
            create_new_DESeq2DataSet(expt_data, baseline_data = Mm_baseline, gender_column = gender_column,
                                     groups = groups, condition_column = condition_column, 
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
                                     groups = groups, condition_column = condition_column, 
                                     match_stages = FALSE, session_obj = session ))

        # experiment data plus stage matched baseline samples
        # design formula includes stage
        groups <- c(groups, 'stage')
        expt_plus_baseline_with_stage_dds <-
          suppressWarnings(
            create_new_DESeq2DataSet(expt_data, baseline_data = Mm_baseline, gender_column = gender_column,
                                     groups = groups, condition_column = condition_column, 
                                     match_stages = TRUE, session_obj = session ))

        return(
          list(
            expt_only_dds = expt_only_dds,
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

  pcs <- reactive({
    if (session$userData[['debug']]) {
      cat("Function: pcs\n")
    }
    pca_info <- pca_info
    if(is.null(pca_info)) {
      return(NULL)
    } else {
      # return names of columns that contribute > 1% of variance
      pca <- pca_info[['pca']]
      return(colnames(pca[['x']])[ pca_info[['propVarPC']] > 0.01 ])
    }
  })
  
  # Change UI options based on input data
  # Components
  observe({
    pc_names <- pcs()
    if (session$userData[['debug']]) {
      cat("Function: UI observer - Components\n")
      print(pc_names)
    }
    if(!is.null(pc_names)) {
      pc_names <- as.list(pc_names)
      names(pc_names) <- pc_names
      # Set the label, choices, and selected item
      updateRadioButtons(session, "x_axis_pc",
                         choices = pc_names,
                         selected = pc_names[[1]]
      )
      updateRadioButtons(session, "y_axis_pc",
                         choices = pc_names,
                         selected = pc_names[[2]]
      )
    }
  })
  
  pca_plot_data <- reactive({
    pca_info <- reactiveValuesToList(pca_info)
    if (is.null(pca_info[['pca']])) {
      return(NULL)
    } else {
      dds_vst <- pca_info[['dds_vst']]
      pca <- pca_info[['pca']]
      plot_data <- cbind( pca[['x']], as.data.frame(colData(dds_vst)) )
      return(plot_data)
    }
  })
  
  # render PCA plot
  output$pca_plot_reduced <- renderPlot({
    plot_data <- pca_plot_data()
    if (is.null(plot_data)) {
      return(NULL)
    } else {
      dds_vst <- pca_info[['dds_vst']]
      if (is.null(dds_vst)) {
        return(NULL)
      } else {
        col_palette <- colour_palette(dds_vst)
        shape_palette <- shape_palette(dds_vst)
        if (session$userData[['debug']]) {
          print(shape_palette)
        }
        pca_plot <- 
          scatterplot_with_fill_and_shape(
            plot_data, input$x_axis_pc, input$y_axis_pc, 
            fill_var = 'stage', fill_palette = col_palette,
            shape_var = 'condition', shape_palette = shape_palette,
            sample_names = input$sample_names)
  
        return(pca_plot)
      }
    }
  })
  
  pca_info_all <- reactiveValues()
  observer_a <- observe({
    expt_data <- exptData()
    deseq_datasets <- deseqDatasets()
    if (is.null(deseq_datasets)) {
      return(NULL)
    } else {
      if (session$userData[['debug']]) {
        print('Function: run_pca_all')
        print('Variance Stabilizing Transform begin...')
      }
      
      expt_plus_all_baseline_dds_vst <- 
        varianceStabilizingTransformation(deseq_datasets[['expt_plus_all_baseline_dds']], blind=TRUE)
      
      if (session$userData[['debug']]) {
        cat('Variance Stabilizing Transform done.\n')
      }
      
      pca <- prcomp( t( assay(expt_plus_all_baseline_dds_vst) ) )
      pca_info_all[['pca']] <- pca
      pca_info_all[['propVarPC']] <- pca$sdev^2 / sum( pca$sdev^2 )
      aload <- abs(pca$rotation)
      pca_info_all[['propVarRegion']] <- sweep(aload, 2, colSums(aload), "/")
      pca_info_all[['dds_vst']] <- expt_plus_all_baseline_dds_vst
    }
  }, priority = -1000)
  
  pca_plot_data_all <- reactive({
    pca_info <- reactiveValuesToList(pca_info_all)
    if (is.null(pca_info[['pca']])) {
      return(NULL)
    } else {
      expt_plus_all_baseline_dds_vst <- pca_info[['dds_vst']]
      pca <- pca_info[['pca']]
      plot_data <- cbind( pca[['x']], as.data.frame(colData(expt_plus_all_baseline_dds_vst)) )
      return(plot_data)
    }
  })
  
  output$pca_plot_all <- renderPlot({
    plot_data <- pca_plot_data_all()
    if (is.null(plot_data)) {
      return(NULL)
    } else {
      expt_plus_all_baseline_dds <- pca_info_all[['dds_vst']]
      if (is.null(expt_plus_all_baseline_dds)) {
        return(NULL)
      } else {
        col_palette <- colour_palette(expt_plus_all_baseline_dds)
        shape_palette <- shape_palette(expt_plus_all_baseline_dds)
        pca_plot <- 
          scatterplot_with_fill_and_shape(
            plot_data, input$x_axis_pc, input$y_axis_pc, 
            fill_var = 'stage', fill_palette = col_palette,
            shape_var = 'condition', shape_palette = shape_palette,
            sample_names = input$sample_names)
        
        return(pca_plot)
      }
    }
  })
  
  # run DESeq2
  deseq_results <- reactive({
    deseq_datasets <- deseqDatasets()
    if (!is.null(deseq_datasets)) {
      deseq_results_3_ways <- 
        overlap_deseq_results( deseq_datasets, 'hom', 'wt', session )
      return(deseq_results_3_ways)
    }
  })
  
  # show results in results tab
  output$results_table <- DT::renderDataTable({
    results_table <- deseq_results()[['merged_results']]
    if (!is.null(results_table)) {
      return(as.data.frame(results_table))
    }
  }, server = TRUE,
  selection = 'single', rownames = FALSE,
  options = list(pageLength = 100))
  
  observeEvent(input$results_table_rows_selected,{
      updateTabsetPanel(session, 'baseline_compare', selected = 'count_plot_panel')
    },
    priority = 1000
  )
  
  output$count_plot_selected_gene <- renderText({
    row_number <- input$results_table_rows_selected
    if (!is.null(row_number)) {
      return(row_number)
    }
  })
}
