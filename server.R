library(shiny)
library(SummarizedExperiment)

# source functions
source(file.path('R', 'load_data.R'))

# load baseline data
# loads object called baseline_data
load(file.path('data', 'Mm_baseline_data.rda'))

# Server logic
server <- function(input, output, session) {
  # set testing and debugging options
  session$userData[['debug']] <- TRUE
  session$userData[['testing']] <- TRUE
  
  # load samples and counts files
  exptData <- reactive({
    if (session$userData[['debug']]) {
      print('Function: mergedCounts')
    }
    if (session$userData[['testing']]) {
      sample_file <- file.path('data', 'Brd2-samples.txt')
      count_file <- file.path('data', 'Brd2-counts.tsv')
      expt_data <-
        load_data(sample_file, count_file, Mm_baseline, session)
    } else{
      # sample_file_info <- input$sample_file
      # count_file_info <- input$count_file
      # if (!is.null(sampleFileInfo) &
      #     !is.null(countFileInfo)) {
      #   exptData <-
      #     load_data(sampleFileInfo$datapath,
      #               countFileInfo$datapath,
      #               input$dataType,
      #               session)
      # } else if (!is.null(dataFileInfo)) {
      #   load(dataFileInfo$datapath)
      # }
      # else{
      #   return(NULL)
      # }
    }
    return(expt_data)
  })
    
  # check columns
  
  # compare to baseline data
  # check sample stages and gene ids
  
  # combine data, getting appropriate baseline samples by stage
  
  # run PCA
  
  # run DESeq2
  
  # overlap the 3 DE lists
  
  # show results in results tab
}
