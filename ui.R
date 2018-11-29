library(shiny)
library(shinycssloaders)
library(shinyjs)
library(shinyBS)

# UI for application
ui <- fluidPage(
  tags$head( includeScript('www/results_table.js') ),
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
                 numericInput(inputId = 'sig_level', label = 'Threshold for statistical significance', 
                              value = 0.05, min = 0, max = 1, step = 0.01, width = NULL),
                 actionButton('analyse_data', 'Analyse Data', icon = NULL, width = NULL),
                 hr(),
                 # options
                 h4('Options'),
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
                 h5("Gender"),
                 checkboxInput("use_gender", 
                               label = 'Use gender information',
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
                              h3('Instructions'),
                              p('Instruction text goes here'))
                   ),
                   fluidRow(
                      h3('Progress'),
                      tags$div(class = "well",
                                bsAlert("progress")),
                      tags$div(class = "well", textOutput('results_text'))
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
                      withSpinner(plotOutput("pca_plot_reduced", height = "480px")),
                      withSpinner(plotOutput("pca_plot_all", height = "480px"))
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
                     ),
                     h4('Downloads'),
                     radioButtons(
                       "plotFormat",
                       label = h5("Plot File"),
                       choices = list('pdf' = 'pdf', 
                                      'eps' = 'eps',
                                      'svg' = 'svg',
                                      'png' = 'png'),
                       selected = 'pdf'
                     ),
                     downloadButton('download_current', 'Download Current Plot'),
                     downloadButton('download_all', 'Download all (pdf)'),
                     downloadButton('download_rda', 'Download rda file of plot')
                   )
                 )
               )
             )
    ),
    tabPanel("Results",
             fluidPage(
               fluidRow(
                 tags$div(class = "well",
                        bsAlert("deseq_progress"))
               ),
               fluidRow(
                 column(width = 3,
                        img(src = 'images/unprocessed_icon_150.png', class = 'results_button', 
                            id = "unprocessed", width = 240),
                        img(src = 'images/mutant_response_icon_150.png', class = 'results_button', 
                            id = "mutant_response", width = 240),
                        img(src = 'images/delay_icon_150.png', class = 'results_button', 
                            id = "delay", width = 240),
                        img(src = 'images/no_delay_icon_150.png', class = 'results_button', 
                            id = "no_delay", width = 240),
                        img(src = 'images/discard_icon_150.png', class = 'results_button', 
                            id = "discard", width = 240),
                        p(class = 'results_button', id = 'all_genes')
                 ),
                 column( width = 9,
                   DT::dataTableOutput(outputId="results_table")
                 )
               )
             )
    ),
    tabPanel("Count Plot", value = 'count_plot_panel',
      plotOutput('count_plot_selected_gene')
    ),
    tabPanel("Help"
    )
  )
)
