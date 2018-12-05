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
  # load baseline data
  # loads object called Mm_baseline
  load(file.path('data', 'Mm_baseline_data.rda'))
  
  # set testing and debugging options
  session$userData[['debug']] <- TRUE
  session$userData[['testing']] <- TRUE
  # session <- list(userData = list(testing = TRUE, debug = TRUE))
  
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
  
  obsConditionLevels <- observe({
    condition_var <- input$condition_var
    expt_data <- isolate(exptData())
    # check levels
    if(valid_condition_column( colData(expt_data)[[condition_var]] )) {
      closeAlert(session, 'invalid_condition_col')
    } else {
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
                          title = 'Input Files', content = conditionMessage(w), style = 'warning')
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
  pca_info <- reactive({
    deseq_datasets <- deseqDatasets()
    if (is.null(deseq_datasets)) {
      return(NULL)
    } else {
      if (session$userData[['debug']]) {
        cat('Function: pca_info\n')
        cat('Variance Stabilizing Transform begin...\n')
      }
      # update progress alert
      closeAlert(session, 'progress_1')
      createAlert(session, anchorId = 'progress', alertId = 'progress_2',
                  content = 'Calculating PCA...', dismiss = FALSE)
      
      # Create a Progress object
      progress <- shiny::Progress$new(session)
      # Make sure it closes when we exit this reactive, even if there's an error
      on.exit(progress$close())
      
      progress$set(message = "Calculating PCA...",
                   detail = 'This will depend on the number of samples', value = 0.2)
      
      dds_vst <- varianceStabilizingTransformation(deseq_datasets[['expt_plus_baseline_dds']], blind=TRUE)
      
      progress$set(value = 0.5)
      
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
  })

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
      if (input$plot_format_reduced == "pdf") {
        pdf(file, paper = "special", height = 7, width = 10) # open the pdf device
      } else if (input$plot_format_reduced == "eps") {
        postscript(file, paper = "special", height = 7, width = 10) # open the postscript device
      } else if (input$plot_format_reduced == "svg") {
        svglite(file, height = 7, width = 10) # open the svg device
      } else if (input$plot_format_reduced == "png") {
        png(file, height = 480, width = 960, res = 100) # open the png device
      } else {
        pdf(file, paper = "special", height = 7, width = 10) # open the pdf device
      }
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
      if (input$plot_format_all == "pdf") {
        pdf(file, paper = "special", height = 7, width = 10) # open the pdf device
      } else if (input$plot_format_all == "eps") {
        postscript(file, paper = "special", height = 7, width = 10) # open the postscript device
      } else if (input$plot_format_all == "svg") {
        svglite(file, height = 7, width = 10) # open the svg device
      } else if (input$plot_format_all == "png") {
        png(file, height = 480, width = 960, res = 100) # open the png device
      } else {
        pdf(file, paper = "special", height = 7, width = 10) # open the pdf device
      }
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
  
  # run DESeq2
  deseq_results <- reactive({
    deseq_datasets <- deseqDatasets()
    if (!is.null(deseq_datasets)) {
      # update progress alert
      createAlert(session, anchorId = 'deseq_progress_1', alertId = 'progress_4',
                  content = 'Running DESeq2...', dismiss = FALSE)
      closeAlert(session, 'deseq_not_started')
      createAlert(session, anchorId = 'deseq_progress_2', alertId = 'progress_5',
                  title = 'DESeq2 Analysis',
                  content = 'Running DESeq2. This may take a while', dismiss = FALSE)
      sig_level <- input$sig_level
      deseq_results_3_ways <- 
        overlap_deseq_results( deseq_datasets, expt_condition = exptCondition(), 
                                ctrl_condition =  ctrlCondition(), sig_level = sig_level,
                                session_obj = session )
      return(deseq_results_3_ways)
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
  
  output$results_table <- DT::renderDataTable({
    results_source <- resultsSource()
    results <- deseq_results()
    if (!is.null(results) & !is.null(results_source)) {
      if (session$userData[['debug']]) {
        cat("Function: results_table\n")
        print(results_source)
        print(class(results))
        print(length(results))
        print(names(results))
        print(class(results[['results_tables']]))
        print(length(results[['results_tables']]))
        print(names(results[['results_tables']]))
      }
      
      results_table <- results[['results_tables']][[results_source]]
      # create datatable, accounting for the different columns
      # in unprocessed vs others
      if (results_source == 'unprocessed') {
        results_dt <- datatable(results_table,
                                selection = 'single', rownames = FALSE,
                                options = list(pageLength = 100)) %>%
          formatStyle(c("padj"), 
                      backgroundColor = styleInterval(0.05, c('white', '#C0C0C0')) ) %>%
          formatStyle(c("log2FC"),
                      valueColumns = c("padj"),
                      backgroundColor = styleInterval(0.05, c('white', '#C0C0C0')) )
        
      } else {
        results_dt <- datatable(results_table,
                                selection = 'single', rownames = FALSE,
                                options = list(pageLength = 100)) %>%
          formatStyle(c("padj.expt_only", "padj.plus_baseline", "padj.with_stage"), 
                      backgroundColor = styleInterval(0.05, c('white', '#C0C0C0')) ) %>%
          formatStyle(c("log2FC.expt_only", "log2FC.plus_baseline", "log2FC.with_stage"),
                      valueColumns = c("padj.expt_only", "padj.plus_baseline", "padj.with_stage"),
                      backgroundColor = styleInterval(0.05, c('white', '#C0C0C0')) )
      }
      return(results_dt)
    }
  }, server = TRUE)
  
  observeEvent(input$results_table_rows_selected,{
      updateTabsetPanel(session, 'baseline_compare', selected = 'count_plot_panel')
    },
    priority = 1000
  )
  
  output$count_plot_selected_gene <- renderPlot({
    row_number <- input$results_table_rows_selected
    results_source <- resultsSource()
    results <- deseq_results()
    if (!is.null(row_number)) {
      closeAlert(session, 'count_plot_instructions')
      # get count data for gene
      results_table <- results[['results_tables']][[results_source]]
      if( session$userData[['debug']] ) {
        print(row_number)
        print(head(results_table))
      }
      gene_id <- results_table[ row_number, 'Gene.ID' ]
      expt_plus_baseline <- results[['plus_baseline_res']]
      counts <- counts(expt_plus_baseline[['deseq']], normalized = TRUE)[ gene_id, ]
      counts_m <- melt(counts, value.name = 'Counts')
      plot_data <- as.data.frame(merge(counts_m, colData(expt_plus_baseline[['deseq']]), by = 'row.names'))
      
      if( session$userData[['debug']] ) {
        print(counts)
        print(head(plot_data))
      }
      
      col_palette <- colour_palette(plot_data$stage)
      shape_palette <- shape_palette(plot_data$condition)
      count_plot <- scatterplot_with_fill_and_shape(plot_data, x_var = 'sample_name', y_var = 'Counts', 
                                                    fill_var = 'stage', fill_palette = col_palette, 
                                                    shape_var = 'condition', shape_palette = shape_palette, 
                                                    sample_names = FALSE) + theme(axis.text.x = element_text(angle = 90))
      
      return(count_plot)
    }
  })
}
