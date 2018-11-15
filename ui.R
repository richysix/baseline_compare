library(shiny)
library(shinycssloaders)

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
                               value = TRUE),
                 hr(),
                 # options
                 # use gender, needs name of gender column as well
                 checkboxInput("use_gender", 
                               label = 'Use gender information',
                               value = TRUE)
               ),
               mainPanel(
                 width = 8,
                 fluidPage(
                   fluidRow(bsAlert("input_file_alert")),
                   fluidRow(withSpinner(textOutput("pca_progress"))),
                   fluidRow(withSpinner(textOutput("deseq_progress")))
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
