library(shiny)
library(shinyjs)
library(shinyBS)
library(DT)
library(SummarizedExperiment)
library(DESeq2)
library(ggplot2)
library(scales)
library(shinyMisc)
library(biovisr)
library(reshape2)
library(svglite)

# source functions
source(file.path('R', 'load_data.R'))
source(file.path('R', 'deseq_functions.R'))
source(file.path('R', 'helper_functions.R'))

# set option to make datatables render NA values as a string
options(htmlwidgets.TOJSON_ARGS = list(na = 'string'))

# set some constants
allowed_conditions <- c('hom', 'het', 'wt', 'mut', 'sib')

# Server logic
server <- function(input, output, session) {
################################################################################
  ## FILE INPUT TAB
  
  # set testing and debugging options
  session$userData[['debug']] <- TRUE
  session$userData[['testing']] <- TRUE
  session$userData[['precomputed']] <- FALSE
  # session <- list(userData = list(testing = TRUE, debug = TRUE))
  
  # find available baseline versions
  get_baseline_files <- function(dir) {
    baseline_files <- grep("Mm_GRCm[0-9]+_e[0-9]+_baseline_data.rda", 
                           list.files(dir), value = TRUE)
    ensembl_names <- gsub("Mm_", "", baseline_files)
    ensembl_names <- gsub("_baseline_data.rda", "", ensembl_names)
    names(baseline_files) <- ensembl_names
    return(baseline_files)
  }
  baselineFiles <- reactive({
    get_baseline_files('data')
  })
  
  # Change UI options based on input data
  # ensembl versions
  observe({
    if (session$userData[['debug']]) {
      cat("Function: UI observer - Ensembl version\n")
    }
    if (!is.null(baselineFiles())) {
      # ensembl version options
      ensembl_versions_options <- as.list(names(baselineFiles()))
      names(ensembl_versions_options) <- names(baselineFiles())
      if ( any(ensembl_versions_options == 'GRCm38_88') ) {
        selected_option <- 'GRCm38_88'
      } else {
        selected_option <- ensembl_versions_options[[1]]
      }
      updateRadioButtons(session, "ensembl_version",
                         choices = ensembl_versions_options,
                         selected = selected_option
      )
    }
  })
  
  # load baseline data
  Mm_baseline <- reactive({
    version_name <- input$ensembl_version
    baseline_files <- baselineFiles()
    if (session$userData[['debug']]) {
      print('Function: Mm_baseline')
      print(version_name)
      print(baseline_files)
    }
    load(file.path('data', baseline_files[version_name]))
    return(eval(as.name(paste0('Mm_', version_name, '_baseline'))))
  })

  exptCondition <- reactiveVal(value = 'mut')
  ctrlCondition <- reactiveVal(value = 'sib')
  
  # create initial alerts
  createAlert(session, anchorId = 'progress', alertId = 'progress_0',
              content = 'Waiting for data upload...', dismiss = FALSE)
  createAlert(session, anchorId = 'count_plot_alert', alertId = 'count_plot_instructions',
              title = 'Instructions', style = 'info',
              content = 'Click on a row in the Results table to see a plot of the normalised counts for each sample')
  createAlert(session, anchorId = 'deseq_progress_2', alertId = 'deseq_not_started',
              title = 'DESeq2 Analysis', style = 'info',
              content = 'DESeq is not running yet. Check the "Files" tab.')
  
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
    # update alerts
    closeAlert(session, 'progress_0')
    createAlert(session, anchorId = 'progress', alertId = 'progress_1',
                content = 'Experiment Data loaded.', dismiss = FALSE)
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
  
  # is the selected condition column valid?
  validConditionColumn <- reactive({
    condition_var <- input$condition_var
    expt_data <- exptData()
    if (!is.null(expt_data)) {
      return( valid_condition_column( colData(expt_data)[[condition_var]] ) )
    }
  })
  
  # create alert if an invalid column is selected
  obsConditionLevels <- observe({
    valid_condition_column <- validConditionColumn()
    if (!is.null(valid_condition_column)) {
      if( valid_condition_column ) {
        closeAlert(session, 'invalid_condition_col')
      } else {
        condition_var <- isolate(input$condition_var)
        expt_data <- isolate(exptData())
        closeAlert(session, 'invalid_condition_col')
        createAlert(
          session, anchorId = 'input_file_alert', dismiss = TRUE,
          alertId = 'invalid_condition_col', title = 'Invalid Condition Column', 
          content = paste('The selected condition column contains', 
                          'entries that not allowed.', 'Valid values are: ',
                          paste0(allowed_conditions, collapse = ', '), '<br>',
                          'Selected column looks like this: ',
                          paste0(colData(expt_data)[[condition_var]], 
                                 collapse = ', ')),
          style = 'danger'
        )
      }
    }
  })
  
  # create DESeq2 data sets
  deseqDatasets <- reactive({
    button_value <- input$analyse_data
    if (button_value > 0) {
      if (session$userData[['precomputed']]) {
        load('data/test_deseq_datasets.rda')
        return(deseq_datasets)
      } else {
        # create alert and stop analysis if an invalid column is selected
        if(!isolate(validConditionColumn())) {
          # close any open alert
          closeAlert(session, 'invalid_condition_col')
          createAlert(
            session, anchorId = 'input_file_alert', dismiss = FALSE,
            alertId = 'invalid_condition_col', title = 'Invalid Condition Column', 
            content = paste('The selected condition column contains', 
                            'entries that not allowed.', 'Valid values are: ',
                            paste0(allowed_conditions, collapse = ', '), '<br>',
                            "DESeq2 analysis can't be started."),
            style = 'danger'
          )
          return(NULL)
        } else {
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
            Mm_baseline <- isolate(Mm_baseline())
            expt_plus_baseline_dds <-
              withCallingHandlers(
                create_new_DESeq2DataSet(expt_data, baseline_data = Mm_baseline, gender_column = gender_column,
                                         groups = groups, condition_column = condition_column, 
                                         match_stages = TRUE, session_obj = session ),
                error = function(e){ stop(e) },
                warning = function(w){
                  print(conditionMessage(w))
                  msg_content <- conditionMessage(w)
                  warning_components <- unlist( strsplit(msg_content, ": ") )
                  alert_title <- warning_components[1]
                  add_button <- FALSE
                  if (alert_title == 'Missing genes in Baseline data') {
                    alert_anchor <- 'input_file_alert_baseline'
                    add_button <- TRUE
                    alert_div <- '#input_file_alert_baseline > .alert'
                    button_name <- 'missing_genes_list_baseline'
                    msg_content <- paste0(warning_components[2],
                                          '. <br>Click the button below to download ',
                                          'a list of the missing gene ids.<br>')
                  } else if (alert_title == 'Missing genes in experimental data') {
                    alert_anchor <- 'input_file_alert_expt'
                    add_button <- TRUE
                    alert_div <- '#input_file_alert_expt > .alert'
                    button_name <- 'missing_genes_list_expt'
                    msg_content <- paste0(warning_components[2],
                                          '. <br>Click the button below to download ',
                                          'a list of the missing gene ids.<br>')
                  } else {
                    alert_anchor <- 'input_file_alert_stages'
                    add_button <- FALSE
                    msg_content <- warning_components[2]
                  }
                  createAlert(session, anchorId = alert_anchor, 
                              title = alert_title, content = msg_content, 
                              style = 'warning')
                  if (add_button) {
                    insertUI(alert_div, where = "beforeEnd",
                             downloadButton(button_name, 'Missing Genes'),
                             multiple = FALSE, immediate = TRUE,
                             session = session)
                  }
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
      }
    }
  })
  
  output$missing_genes_list_expt <- downloadHandler(
    filename = function() {
      paste(Sys.Date(), 'missing_genes-expt.tsv', sep = '.')
    },
    content = function(file) {
      expt_data <- isolate(exptData())
      baseline_only <- setdiff(rownames(Mm_baseline()), rownames(expt_data))
      
      write.table(
        baseline_only, file = file, quote = FALSE,
        col.names = FALSE, row.names = FALSE, sep = "\t"
      )
    },
    contentType = 'text/tsv'
  )
  
  output$missing_genes_list_baseline <- downloadHandler(
    filename = function() {
      paste(Sys.Date(), 'missing_genes-baseline.tsv', sep = '.')
    },
    content = function(file) {
      expt_data <- isolate(exptData())
      expt_only_genes <- setdiff(rownames(expt_data), rownames(Mm_baseline()))
      
      write.table(
        expt_only_genes, file = file, quote = FALSE,
        col.names = FALSE, row.names = FALSE, sep = "\t"
      )
    },
    contentType = 'text/tsv'
  )
  
  # run PCA with just match baseline samples for speed
  pca_info <- reactive({
    deseq_datasets <- deseqDatasets()
    if (is.null(deseq_datasets)) {
      return(NULL)
    } else {
      if (session$userData[['precomputed']]) {
        load('data/test_pca_info.rda')
        return(info)
      } else {
        if (session$userData[['debug']]) {
          cat('Function: pca_info\n')
          cat('Variance Stabilizing Transform begin...\n')
        }
        # update progress alert
        closeAlert(session, 'progress_1')
        closeAlert(session, 'progress_3')
        createAlert(session, anchorId = 'progress', alertId = 'progress_2',
                    content = 'Calculating PCA...', dismiss = FALSE)
        
        # Create a Progress object
        progress <- shiny::Progress$new(session)
        # Make sure it closes when we exit this reactive, even if there's an error
        on.exit(progress$close())
        
        progress$set(message = "Calculating PCA...",
                     detail = 'This will depend on the number of samples', value = 0.25)
        
        dds_vst <- varianceStabilizingTransformation(deseq_datasets[['expt_plus_baseline_dds']], blind=TRUE)
        
        progress$set(value = 0.4)
        
        if (session$userData[['debug']]) {
          cat('Variance Stabilizing Transform done.\n')
        }
        
        info <- list()
        pca <- prcomp( t( assay(dds_vst) ) )
        info[['subset']][['pca']] <- pca
        info[['subset']][['propVarPC']] <- pca$sdev^2 / sum( pca$sdev^2 )
        aload <- abs(pca$rotation)
        info[['subset']][['propVarRegion']] <- sweep(aload, 2, colSums(aload), "/")
        info[['subset']][['dds_vst']] <- dds_vst
        
        if (session$userData[['debug']]) {
          cat('Function: pca_info, all samples\n')
          cat('Variance Stabilizing Transform begin...\n')
        }
        
        expt_plus_all_baseline_dds_vst <- 
          varianceStabilizingTransformation(deseq_datasets[['expt_plus_all_baseline_dds']], blind=TRUE)
        
        if (session$userData[['debug']]) {
          cat('Variance Stabilizing Transform done.\n')
        }
        
        pca <- prcomp( t( assay(expt_plus_all_baseline_dds_vst) ) )
        info[['all']][['pca']] <- pca
        info[['all']][['propVarPC']] <- pca$sdev^2 / sum( pca$sdev^2 )
        aload <- abs(pca$rotation)
        info[['all']][['propVarRegion']] <- sweep(aload, 2, colSums(aload), "/")
        info[['all']][['dds_vst']] <- expt_plus_all_baseline_dds_vst
        
        progress$set(value = 1)
        
        # update progress alert
        closeAlert(session, 'progress_2')
        createAlert(session, anchorId = 'progress', alertId = 'progress_3',
                    content = 'PCA Finished', dismiss = FALSE)
        
        return(info)
      }
    }
  })

################################################################################
  ## PCA TAB
  pcs <- reactive({
    info <- pca_info()
    if (session$userData[['debug']]) {
      cat("Function: pcs\n")
    }
    if(is.null(info)) {
      return(NULL)
    } else {
      # return names of columns that contribute > 1% of variance
      pca <- info[['subset']][['pca']]
      return(colnames(pca[['x']])[ info[['subset']][['propVarPC']] > 0.01 ])
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
  
  pcaPlotDataReduced <- reactive({
    pca_info <- pca_info()
    if (is.null(pca_info)) {
      return(NULL)
    } else {
      dds_vst <- pca_info[['subset']][['dds_vst']]
      pca <- pca_info[['subset']][['pca']]
      plot_data <- cbind( pca[['x']], as.data.frame(colData(dds_vst)) )
      return(plot_data)
    }
  })
  
  pcaPlotReduced <- reactive({
    plot_data <- pcaPlotDataReduced()
    if(is.null(plot_data)) {
      return(NULL)
    } else {
      col_palette <- colour_palette(plot_data[['stage']])
      shape_palette <- shape_palette(plot_data[['condition']])
      if (session$userData[['debug']]) {
        print(col_palette)
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
  })
  
  # render PCA plot
  output$pca_plot_reduced <- renderPlot({
    return(pcaPlotReduced())
  })
  
  # for downloading the plot as a pdf/png
  output$download_current_pca_reduced <- downloadHandler(
    filename = function() {
      paste('pca_plot', Sys.Date(), input$plot_format_reduced, sep = '.')
    },
    content = function(file) {
      open_graphics_device(input$plot_format_reduced, file)
      print(pcaPlotReduced())
      dev.off()  # close device
    },
    contentType = paste0('image/', input$plot_format_reduced)
  )
  
  # download an rda file of the current plot
  output$download_rda_pca_reduced <- downloadHandler(
    filename = function() {
      paste('pca_plot', Sys.Date(), 'rda', sep = '.')
    },
    content = function(file) {
      pca_plot <- pcaPlotReduced()
      save(pca_plot, file = file)
    }
  )
  
  # for downloading a pdf of each component plotted against the next one
  output$download_all_pca_reduced <- downloadHandler(
    filename = function() {
      paste('pca_plot-allPCs', Sys.Date(), 'pdf', sep = '.')
    },
    content = function(file) {
      pdf(file, paper = "special", height = 7, width = 10) # open the pdf device
      
      plot_data <- pcaPlotDataReduced()
      col_palette <- colour_palette(plot_data[['stage']])
      shape_palette <- shape_palette(plot_data[['condition']])
      components <- pcs()
      for (i in seq.int(length(components) - 1)) {
        first <- components[i]
        second <- components[i + 1]
        pca_plot <- 
          scatterplot_with_fill_and_shape(
            plot_data, first, second, 
            fill_var = 'stage', fill_palette = col_palette,
            shape_var = 'condition', shape_palette = shape_palette,
            sample_names = input$sample_names)
        print(pca_plot)
      }
      dev.off()  # close device
    },
    contentType = 'image/pdf'
  )
  
  pcaPlotDataAll <- reactive({
    pca_info <- pca_info()
    if (is.null(pca_info)) {
      return(NULL)
    } else {
      expt_plus_all_baseline_dds <- pca_info[['all']][['dds_vst']]
      pca <- pca_info[['all']][['pca']]
      plot_data <- cbind( pca[['x']], as.data.frame(colData(expt_plus_all_baseline_dds)) )
      return(plot_data)
    }
  })
  
  pcaPlotAll <- reactive({
    plot_data <- pcaPlotDataAll()
    if (is.null(plot_data)) {
      return(NULL)
    } else {
      col_palette <- colour_palette(plot_data[['stage']])
      shape_palette <- shape_palette(plot_data[['condition']])
      pca_plot <- 
        scatterplot_with_fill_and_shape(
          plot_data, input$x_axis_pc, input$y_axis_pc, 
          fill_var = 'stage', fill_palette = col_palette,
          shape_var = 'condition', shape_palette = shape_palette,
          sample_names = input$sample_names)
      return(pca_plot)
    }
  })
  
  output$pca_plot_all <- renderPlot({
    return(pcaPlotAll())
  })
  
  # for downloading the plot as a pdf/png
  output$download_current_pca_all <- downloadHandler(
    filename = function() {
      paste('pca_plot_all', Sys.Date(), input$plot_format_all, sep = '.')
    },
    content = function(file) {
      open_graphics_device(input$plot_format_all, file)
      print(pcaPlotAll())
      dev.off()  # close device
    },
    contentType = paste0('image/', input$plot_format_all)
  )
  
  # download an rda file of the current plot
  output$download_rda_pca_all <- downloadHandler(
    filename = function() {
      paste('pca_plot_all', Sys.Date(), 'rda', sep = '.')
    },
    content = function(file) {
      pca_plot <- pcaPlotAll()
      save(pca_plot, file = file)
    }
  )
  
  # for downloading a pdf of each component plotted against the next one
  output$download_all_pca_all <- downloadHandler(
    filename = function() {
      paste('pca_plot_all-allPCs', Sys.Date(), 'pdf', sep = '.')
    },
    content = function(file) {
      pdf(file, paper = "special", height = 7, width = 10) # open the pdf device
      
      plot_data <- pcaPlotDataAll()
      col_palette <- colour_palette(plot_data[['stage']])
      shape_palette <- shape_palette(plot_data[['condition']])
      components <- pcs()
      for (i in seq.int(length(components) - 1)) {
        first <- components[i]
        second <- components[i + 1]
        pca_plot <- 
          scatterplot_with_fill_and_shape(
            plot_data, first, second, 
            fill_var = 'stage', fill_palette = col_palette,
            shape_var = 'condition', shape_palette = shape_palette,
            sample_names = input$sample_names)
        print(pca_plot)
      }
      dev.off()  # close device
    },
    contentType = 'image/pdf'
  )
  
################################################################################
  ## RESULTS TAB
  # run DESeq2
  deseq_results <- reactive({
    deseq_datasets <- deseqDatasets()
    if (!is.null(deseq_datasets)) {
      if (session$userData[['precomputed']]) {
        load('data/test_deseq_results_object.rda')
        return(deseq_results_3_ways)
      } else {
        # update progress alert
        closeAlert(session, 'progress_6')
        closeAlert(session, 'progress_7')
        createAlert(session, anchorId = 'deseq_progress_1', alertId = 'progress_4',
                    content = 'Running DESeq2...', dismiss = FALSE)
        closeAlert(session, 'deseq_not_started')
        createAlert(session, anchorId = 'deseq_progress_2', alertId = 'progress_5',
                    title = 'DESeq2 Analysis',
                    content = 'Running DESeq2. This may take a while', dismiss = FALSE)
        sig_level <- isolate(input$sig_level)
        shinyjs::removeCssClass(id = 'deseq_output', class = "hidden")
        # empty console output
        shinyjs::html("deseq_console_ouput", html = '')
        
        # Create a Progress object
        progress <- shiny::Progress$new(session)
        # Make sure it closes when we exit this reactive, even if there's an error
        on.exit(progress$close())
        
        progress$set(message = "Running DESeq2...", value = 0.05)
        session$userData[['progress_object']] <- progress
        
        deseq_results_3_ways <- 
          withCallingHandlers(
            overlap_deseq_results( deseq_datasets, expt_condition = exptCondition(), 
                                  ctrl_condition =  ctrlCondition(), sig_level = sig_level,
                                  session_obj = session ),
            message = function(m) {
              message_text <- m$message
              if (!grepl("Experimental Data", message_text) ){ 
                message_text <- paste0("  ", message_text)
              }
              if (grepl("replacing outliers", message_text)) {
                message_text <- gsub("\n", "\n  ", message_text)
                message_text <- gsub("  $", "", message_text)
              }
              shinyjs::html("deseq_console_ouput", html = message_text, add = TRUE)
            }
          )
        progress$set(value = 1)
        return(deseq_results_3_ways)
      }
    }
  })
  
  output$results_text <- renderText({
    results <- deseq_results()
    if (!is.null(results)) {
      # save(results, file = 'data/deseq_results_object.RData')
      closeAlert(session, 'progress_4')
      closeAlert(session, 'progress_5')
      createAlert(session, anchorId = 'deseq_progress_1', alertId = 'progress_6',
                  content = 'DESeq2 Finished', dismiss = FALSE)
      createAlert(session, anchorId = 'deseq_progress_2', alertId = 'progress_7',
                  title = 'DESeq Analysis',
                  content = 'DESeq2 Finished', dismiss = TRUE)
      return('')
    }
  })
  
  # show results in results tab
  resultsSource <- reactiveVal()
  observeEvent(input$mutant_response, { resultsSource('mutant_response') })
  observeEvent(input$delay, { resultsSource('delay') })
  observeEvent(input$no_delay, { resultsSource('no_delay') })
  observeEvent(input$discard, { resultsSource('discard') })
  observeEvent(input$unprocessed, { resultsSource('unprocessed') })
  observeEvent(input$all_genes, { resultsSource('all_genes') })
  observeEvent(resultsSource(), {
    session$sendCustomMessage("selected_results_button", resultsSource())
  })
    
  resultsTable <- reactive({
    results_source <- resultsSource()
    results <- deseq_results()
    if (!is.null(results) & !is.null(results_source)) {
      if (session$userData[['debug']]) {
        cat("Function: results_table\n")
        cat( sprintf('Results: Source = %s, Class = %s, Length = %s\n', 
                    results_source, class(results), length(results)) )
        print(names(results))
        cat( sprintf('Results Tables: Class = %s, Length = %s\n', 
               class(results[['results_tables']]), 
               length(results[['results_tables']])) )
        print(names(results[['results_tables']]))
      }
      
      return( results[['results_tables']][[results_source]] )
    } else {
      return(NULL)
    }
  })
  
  # number of significant genes
  output$num_sig_genes <- renderText({
    results_source <- resultsSource()
    results_table <- resultsTable()
    if (!is.null(results_table)) {
      if (results_source == 'all_genes') {
        return(NULL)
      } else {
        text <- paste0(nrow(results_table), ' significant genes')
        return(text)
      }
    }
  })
  
  # custom header for datatables
  table_header_3ways <- htmltools::withTags(table(
    class = 'display',
    thead(
      tr(
        th(rowspan = 2, 'Gene.ID'),
        th(rowspan = 2, 'Name'),
        th(colspan = 2, 'Experiment Data Only'),
        th(colspan = 2, 'Plus Baseline'),
        th(colspan = 2, 'Plus Baseline with Stage'),
        th(rowspan = 2, 'Results Set'),
        th(rowspan = 2, 'Chr'),
        th(rowspan = 2, 'Start'),
        th(rowspan = 2, 'End'),
        th(rowspan = 2, 'Strand')
      ),
      tr(
        lapply(rep(c('log2FC', 'padj'), 3), th)
      )
    )
  ))
  
  output$results_table <- DT::renderDataTable({
    results_table = resultsTable()
    results_source <- resultsSource()
    if(!is.null(results_table)) {
      # round log2fc and padj columns
      col_names <- names(results_table)
      results_table <- 
        as.data.frame(
          lapply(names(results_table),
              function(col_name){
                if (grepl("log2FC", col_name)) {
                  return(as.numeric(sprintf('%.3f', results_table[[col_name]])))
                } else if (grepl("padj", col_name)) {
                  return(as.numeric(sprintf('%.3g', results_table[[col_name]])))
                } else {
                  return(results_table[[col_name]])
                }
              }
            )
        )
      names(results_table) <- col_names
      # create datatable, accounting for the different columns
      # in unprocessed vs others
      data_table_options <- list(
        pageLength = 100,
        order = list(list(3, 'asc')),
        rowCallback = JS("function( row, data, dataIndex ) {",
                          "colourCells( row, data, dataIndex )}") #colourCells is defined in results_table.js
      )
      if (results_source == 'unprocessed') {
        results_dt <- 
          datatable(results_table, 
            colnames = c("Gene.ID", "Name", "log2FC", "padj", "Results Set", 
                         "Chr", "Start", "End", "Strand" ),
            selection = 'single', rownames = FALSE,
            options = data_table_options ) %>%
          formatStyle(c("padj.expt_only"), 
                      backgroundColor = styleInterval(0.05, c('white', '#C0C0C0')) ) %>%
          formatStyle(c("log2FC.expt_only"),
                      valueColumns = c("padj.expt_only"),
                      backgroundColor = styleInterval(0.05, c('white', '#C0C0C0')) )
        
      } else {
        results_dt <- 
          datatable(results_table, container = table_header_3ways,
            selection = 'single', rownames = FALSE,
            options = data_table_options ) %>%
          formatStyle(c("padj.expt_only", "padj.plus_baseline", "padj.with_stage"), 
                      backgroundColor = styleInterval(0.05, c('white', '#C0C0C0')) ) %>%
          formatStyle(c("log2FC.expt_only", "log2FC.plus_baseline", "log2FC.with_stage"),
                      valueColumns = c("padj.expt_only", "padj.plus_baseline", "padj.with_stage"),
                      backgroundColor = styleInterval(0.05, c('white', '#C0C0C0')) )
      }
      return(results_dt)
    }
  }, server = TRUE)
  
  # for downloading a file of the results for all genes
  output$download_results_all <- downloadHandler(
    filename = function() {
      paste(Sys.Date(), paste(exptCondition(), ctrlCondition(), sep = '_vs_'), 
            'all.tsv', sep = '.')
    },
    content = function(file) {
      deseq_results <- deseq_results()
      results_table <- deseq_results[['results_tables']][['all_genes']]
      write.table(
        results_table, file = file, quote = FALSE,
        col.names = TRUE, row.names = FALSE, sep = "\t"
      )
    },
    contentType = 'text/tsv'
  )
  
  output$download_results_rda <- downloadHandler(
    filename = function() {
      paste('deseq_results', 
            paste(exptCondition(), ctrlCondition(), sep = '_vs_'),
            Sys.Date(), 'rda', sep = '.')
    },
    content = function(file) {
      deseq_results <- deseq_results()
      save(deseq_results, file = file)
    }
  )
  
  observeEvent(input$results_table_rows_selected,{
      updateTabsetPanel(session, 'baseline_compare', 
                        selected = 'count_plot_panel')
    },
    priority = 1000
  )
  
################################################################################
  ## COUNT PLOT TAB
  countPlot <- reactive({
    row_number <- input$results_table_rows_selected
    if (!is.null(row_number)) {
      closeAlert(session, 'count_plot_instructions')
      
      results <- deseq_results()
      results_table <- resultsTable()
      # get gene id by row number
      gene_id <- results_table[ row_number, 'Gene.ID' ]
      gene_name <- results_table[ row_number, 'Name' ]
      # get counts for gene id
      expt_plus_baseline <- results[['plus_baseline_res']]
      counts <- 
        counts(expt_plus_baseline[['deseq']], normalized = TRUE)[ gene_id, ]
      counts_m <- melt(counts, value.name = 'Counts')
      plot_data <- as.data.frame(merge(counts_m, 
                                       colData(expt_plus_baseline[['deseq']]), 
                                       by = 'row.names'))
      
      if( session$userData[['debug']] ) {
        print(row_number)
        print(head(results_table))
        print(counts)
        print(head(plot_data))
      }
      
      col_palette <- colour_palette(plot_data$stage)
      shape_palette <- shape_palette(plot_data$condition)
      count_plot <- 
        scatterplot_with_fill_and_shape(
          plot_data, x_var = 'sample_name', y_var = 'Counts', 
          fill_var = 'stage', fill_palette = col_palette, 
          shape_var = 'condition', shape_palette = shape_palette, 
          sample_names = FALSE
        ) +
        labs(x = 'Sample', y = 'Normalised Counts',
             fill = 'Stage', shape = 'Condition',
             title = paste0(gene_name, ' (', gene_id, ')') ) +
        theme(axis.text.x = element_text(angle = 90),
              plot.title = element_text(face = 'bold', hjust = 0.5))
      
      return(count_plot)
    }
  })
  
  output$count_plot_selected_gene <- renderPlot({
    return(countPlot())
  })
  
  # for downloading the plot as a pdf/png
  output$download_current_count_plot <- downloadHandler(
    filename = function() {
      paste(Sys.Date(), paste(exptCondition(), ctrlCondition(), sep = '_vs_'), 
            'counts', input$plot_format_counts, sep = '.')
    },
    content = function(file) {
      open_graphics_device(input$plot_format_counts, file)
      print(countPlot())
      dev.off()  # close device
    },
    contentType = paste0('image/', input$plot_format_counts)
  )
}
