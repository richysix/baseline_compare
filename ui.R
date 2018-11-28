library(shiny)
library(shinycssloaders)
library(shinyjs)
library(shinyBS)

# UI for application
ui <- fluidPage(
  navbarPage(
    "Baseline CompaRe", id = 'baseline_compare',
    selected = "file_input",
    tabPanel("Files", value = "file_input",
             sidebarLayout(
               sidebarPanel(
                 width = 4,
                 fileInput('sample_file', 'Load Sample File'),
                 fileInput('count_file', 'Load Count File'),
                 checkboxInput("demo_data", 
                               label = 'Use Demo data',
                               value = FALSE),
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
                   fluidRow(textOutput("pca_progress")),
                   fluidRow(textOutput("deseq_progress")),
                   fluidRow(bsAlert("input_file_alert"))
                 )
               )
            )
    ),
    tabPanel("PCA", value = "pca_panel",
             fluidPage(
               fluidRow(
                 sidebarLayout(
                   mainPanel(
                     fluidPage(
                       withSpinner(plotOutput("pca_plot_reduced", height = "480px")),
                       withSpinner(plotOutput("pca_plot_all", height = "480px"))
                     ),
                     width = 9
                   ),
                   sidebarPanel(
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
                     downloadButton('download_rda', 'Download rda file of plot'),
                     width = 3
                   )
                 )
               )
             )
    ),
    tabPanel("Results",
             DT::dataTableOutput(outputId="results_table")
    ),
    tabPanel("Count Plot", value = 'count_plot_panel',
      plotOutput('count_plot_selected_gene')
    ),
    tabPanel("Help"
    )
  )
)
