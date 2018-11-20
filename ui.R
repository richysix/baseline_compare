library(shiny)
library(shinycssloaders)
library(shinyjs)
library(shinyBS)

# UI for application
ui <- fluidPage(
  navbarPage(
    "Baseline CompaRe",
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
                   fluidRow(bsAlert("input_file_alert")),
                   fluidRow(textOutput("pca_progress")),
                   fluidRow(textOutput("deseq_progress"))
                 )
               )
            )
    ),
    tabPanel("PCA", value = "pca_panel",
             withSpinner(plotOutput("pca_plot"))
    ),
    tabPanel("Results"
    ),
    tabPanel("Count Plot"
    ),
    tabPanel("Help"
    )
  )
)
