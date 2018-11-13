library(shiny)
library(SummarizedExperiment)

# load baseline data
# loads object called baseline_data
load(file.path('data', 'Mm_baseline_data.rda'))

# Server logic
server <- function(input, output, session) {
  # set testing and debugging options
  session$userData[['debug']] <- TRUE
  session$userData[['testing']] <- TRUE
  
  # load samples and counts files
  # check columns
  
  # compare to baseline data
  # check sample stages and gene ids
  
  # combine data, getting appropriate baseline samples by stage
  
  # run PCA
  
  # run DESeq2
  
  # overlap the 3 DE lists
  
  # show results in results tab
}
