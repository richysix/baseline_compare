library(shiny)

# UI for application
ui <- fluidPage(
  navbarPage(
    "Baseline CompaRe",
    tabPanel("Files"
    ),
    tabPanel("PCA",
             plotOutput("pca_plot")
    ),
    tabPanel("Results"
    ),
    tabPanel("Count Plot"
    ),
    tabPanel("Help"
    )
  )
)
