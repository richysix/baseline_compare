library(shiny)
library(shinycssloaders)
library(shinyjs)
library(shinyBS)

# UI for application
ui <- fluidPage(
  useShinyjs(),
  tags$head( includeScript('www/results_table.js'),
             includeScript('www/files.js'),
             includeCSS('www/baseline_compare.css') ),
  navbarPage(
    "Baseline CompaRe", id = 'baseline_compare',
    selected = "file_input",
    tabPanel("Files", value = "file_input",
             sidebarLayout(
               position = 'right',
               sidebarPanel(
                 width = 4,
                 fileInput('sample_file', 'Load Sample File'),
                 fileInput('count_file', 'Load Count File'),
                 checkboxInput("demo_data", 
                               label = 'Use Demo data',
                               value = FALSE),
                 actionButton('analyse_data', 'Analyse Data'),
                 p('This button will rerun the analysis every time it is clicked'),
                 hr(),
                 # options
                 h4('Options'),
                 numericInput(inputId = 'sig_level', label = 'Threshold for statistical significance', 
                              value = 0.05, min = 0, max = 1, step = 0.01, width = NULL),
                 radioButtons(
                   "condition_var",
                   label = h5("Condition"),
                   choices = list(
                     "None" = "None",
                     "condition" = 'condition'
                   ),
                   selected = 'condition'
                 ),
                 # use gender, needs name of gender column as well
                 br(),
                 h5("Sex"),
                 checkboxInput("use_gender", 
                               label = 'Use sex information in DESeq',
                               value = TRUE),
                 radioButtons(
                   "sex_var",
                   label = h6("Sex Variable Name"),
                   choices = list(
                     "None" = "None",
                     "sex" = 'sex'
                   ),
                   selected = 'sex'
                 )
               ),
               mainPanel(
                 width = 8,
                 fluidPage(
                   fluidRow(
                     tags$div(class = "well",
                        includeHTML('www/front_page_intro.html')
                     )
                   ),
                   fluidRow(
                      h3('Progress'),
                      tags$div(class = "well",
                        bsAlert("progress"),
                        bsAlert("deseq_progress_1"),
                        div(id = "deseq_output", class = "hidden",
                            h6('DESeq2 output'),
                            pre(id = "deseq_console_ouput")
                        ),
                        tags$div(id = 'deseq_results_text',
                          textOutput('results_text')
                        )
                      )
                   ),
                   fluidRow(
                      bsAlert("input_file_alert")
                   )
                 )
               )
            )
    ),
    tabPanel("PCA", value = "pca_panel",
             fluidPage(
               fluidRow(h2('Principal Component Analysis')),
               fluidRow(
                 sidebarLayout(
                   mainPanel(
                      width = 9,
                      h3('PCA: Baseline sample for experimental stages included'),
                      withSpinner(plotOutput("pca_plot_reduced", height = "480px")),
                      h4('Download'),
                      fluidRow(
                        column( width = 4,
                          radioButtons(
                          "plot_format_reduced",
                          label = NULL,
                          choices = list('pdf' = 'pdf', 
                                         'eps' = 'eps',
                                         'svg' = 'svg',
                                         'png' = 'png'),
                          selected = 'pdf',
                          inline = TRUE
                          ),
                          downloadButton('download_current_pca_reduced', 'Download Current Plot')
                        ),
                        column( width = 8,
                          downloadButton('download_rda_pca_reduced', 'Download rda file of plot'),
                          downloadButton('download_all_pca_reduced', 'Download all PCs (pdf)')
                        )
                      ),
                      h3('PCA: All Baseline samples included'),
                      withSpinner(plotOutput("pca_plot_all", height = "480px")),
                      h4('Download'),
                      fluidRow(
                        column( width = 4,
                                radioButtons(
                                  "plot_format_all",
                                  label = NULL,
                                  choices = list('pdf' = 'pdf', 
                                                 'eps' = 'eps',
                                                 'svg' = 'svg',
                                                 'png' = 'png'),
                                  selected = 'pdf',
                                  inline = TRUE
                                ),
                                downloadButton('download_current_pca_all', 'Download Current Plot')
                        ),
                        column( width = 8,
                                downloadButton('download_rda_pca_all', 'Download rda file of plot'),
                                downloadButton('download_all_pca_all', 'Download all PCs (pdf)')
                        )
                      ),
                      hr()
                   ),
                   sidebarPanel(
                     width = 3,
                     h5("PCA Options"),
                     # checkbox for displaying sample names
                     checkboxInput("sample_names", 
                                   label = 'Sample Names',
                                   value = FALSE),
                     h5("Axes"),
                     # x axis buttons
                     radioButtons(
                       "x_axis_pc",
                       label = h6("X axis component"),
                       choices = list(
                         "PC1" = "PC1",
                         "PC2" = "PC2"
                       ),
                       selected = 1
                     ),
                     # y axis buttons
                     radioButtons(
                       "y_axis_pc",
                       label = h6("Y axis component"),
                       choices = list(
                         "PC1" = "PC1",
                         "PC2" = "PC2"
                       ),
                       selected = 2
                     )
                   )
                 )
               )
             )
    ),
    tabPanel("Results",
             fluidPage(
               fluidRow(bsAlert("deseq_progress_2")),
               fluidRow(
                 column(width = 3,
                        actionButton('mutant_response', 'Mutant Response', 
                                     width = '200px', class = "results-button"),
                        actionButton('delay', 'Delay', width = '200px', 
                                     class = "results-button"),
                        actionButton('no_delay', 'No Delay', width = '200px', 
                                     class = "results-button"),
                        actionButton('discard', 'Discard', width = '200px', 
                                     class = "results-button"),
                        hr(),
                        actionButton('unprocessed', 'Experiment Samples only', 
                                     width = '200px', class = "results-button"),
                        actionButton('all_genes', 'All Genes', width = '200px', 
                                     class = "results-button"),
                        h4('Download'),
                        downloadButton('download_results_all',
                                       'Download Results (All Genes)'),
                        downloadButton('download_results_rda',
                                       'Download Results Objects (.rda)')
                 ),
                 column( width = 9,
                  textOutput('num_sig_genes'),  
                  DT::dataTableOutput(outputId="results_table")
                 )
               )
             )
    ),
    tabPanel("Count Plot", value = 'count_plot_panel',
             bsAlert("count_plot_alert"),
             sidebarLayout(
               position = 'right',
               sidebarPanel(
                 width = 3,
                 h4('Download'),
                 radioButtons(
                   "plot_format_counts",
                   label = 'Plot Format',
                   choices = list('pdf' = 'pdf', 
                                  'eps' = 'eps',
                                  'svg' = 'svg',
                                  'png' = 'png'),
                   selected = 'pdf'
                 ),
                 downloadButton('download_current_count_plot', 'Download Current Plot')
               ),
               mainPanel(
                 width = 9,
                 plotOutput('count_plot_selected_gene')
               )
             )
    ),
    tabPanel("Help", value = 'help_panel',
             includeMarkdown('README.md')
    )
  )
)
