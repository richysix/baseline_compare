library(shiny)
library(shinycssloaders)

# UI for application
ui <- fluidPage(
  navbarPage(
    "Baseline CompaRe",
    tabPanel("Files"
    ),
    tabPanel("PCA",
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
