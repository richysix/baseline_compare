library(shiny)
library(shinycssloaders)

# UI for application
ui <- fluidPage(
  navbarPage(
    "Baseline CompaRe",
    selected = "file_input",
    tabPanel("Files", value = "file_input",
             fileInput('sample_file', 'Load Sample File'),
             fileInput('count_file', 'Load Count File'),
             checkboxInput("test_data", 
                           label = 'Use test data',
                           value = FALSE),
             hr()
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
