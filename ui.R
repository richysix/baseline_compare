library(shiny)
library(shinycssloaders)

# UI for application
ui <- fluidPage(
  navbarPage(
    "Baseline CompaRe",
    selected = "pca_panel",
    tabPanel("Files"
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
