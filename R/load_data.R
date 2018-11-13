#' load_data
#'
#' \code{load_data} returns a Summarized Experiment object produced from the supplied sample and count files 
#'
#'    The first column of the sample file must match the sample names in the column names of the count file.
#'    A warning will be produced detailing sample names present in the sample file but not the count file.
#'    For this, the sample file must also have a column labelled "condition" for the DESeq2 design formula.
#'    
#'    If the samples file contains a column labelled "sampleName" the samples will be renamed.
#'    
#' @param sample_file character - name of the sample file
#' @param count_file character - name of the count file
#' @param baseline_data baseline data object
#' @param session session_object
#' 
#' @return matrix
#'
#' @examples
#' load_data( 'samples.txt', 'all.tsv', Mm_e9-5_baseline, session_obj )
#'
#' @export
#'
load_data <- function( sample_file, count_file, baseline_data, session ){
  # read in sample info
  sample_info <- load_sample_data(sample_file, session)
  ## TO DO
  # make sure stage column is compatible with baseline data
  
  # load count data
  expt_count_data <- load_count_data(count_file, session)
  ## TO DO
  # check for missing genes and warn

  # make a list with all the data in
  all_data <- list(
    counts = as.matrix(countData),
    rowData = data.frame(data[, !grepl('count', names(data))],
                         row.names = rownames(countData)),
    colData = samples
  )
  return(all_data)
}

#' load_sample_data
#'
#' \code{load_sample_data} Reads in the sample data and retuns it
#'
#'    The first column is used as row names.
#'    A sample_name column is created using the row names so that 
#'    the sample data can be joined to the PCA data.
#'    The levels of the column variables are set to the order that 
#'    they appear in the file.
#'    
#' @param sample_file   character - name of the sample file
#' @param session session_object
#' 
#' @return data.frame
#'
#' @examples
#' load_sample_data( 'samples.txt', session_obj )
#'
#' @export
#'
load_sample_data <- function(sample_file, session){
  sample_info <- read.delim(sample_file, row.names = 1,
                            check.names = FALSE)
  sample_info$sample_name <- rownames(sample_info)
  for (column_name in colnames(sample_info)) {
    if (class(sample_info[[column_name]]) != 'factor') {
      next
    }
    # set levels to the order in which they appear in the file
    sample_info[[column_name]] <-
      factor(sample_info[[column_name]],
             levels = unique(sample_info[[column_name]]))
    # unless levels are 'wt', 'het' and 'hom'
    if (all( Reduce('|', 
                    lapply(c('wt', 'het', 'hom'), 
                           function(gt){ levels(sample_info[[column_name]]) == gt }) ) ) ) {
      sample_info[[column_name]] <-
        factor(sample_info[[column_name]],
               levels = c('wt', 'het', 'hom'))
    }
  }
  
  if (session$userData[['debug']]) {
    print(head(sample_info))
  }
  
  return(sample_info)
}

#' load_count_data
#'
#' \code{load_count_data} Reads in the sample data and retuns it
#'
#'    The first column is used as row names.
#'    A sample_name column is created using the row names so that 
#'    the sample data can be joined to the PCA data.
#'    The levels of the column variables are set to the order that 
#'    they appear in the file.
#'    
#' @param count_file   character - name of the sample file
#' @param session session_object
#' 
#' @return data.frame
#'
#' @examples
#' load_count_data( 'counts.txt', session_obj )
#'
#' @export
#'
load_count_data <- function(count_file, session){
  # Read in count data
  data <- read.delim(count_file, check.names=FALSE)
  
  # Support different column names
  names(data)[names(data) == 'Gene ID']      <- 'Gene.ID'
  names(data)[names(data) == 'ID']      <- 'Gene.ID'
  names(data)[ grepl("e[0-9][0-9] Ensembl Gene ID", names(data), perl=TRUE ) ]    <- 'Gene.ID'
  
  ## TO DO
  # check required columns exist
  
  # convert gene id to character
  if (any(grepl('Gene.ID', names(data)))) {
    data[['Gene.ID']] <- as.character(data[['Gene.ID']])
  }
  rownames(data) <- data[ , "Gene.ID" ]
  
  return(data)
}

#' merge_with_baseline
#'
#' \code{merge_with_baseline} Merge the expt and baseline data
#'
#'    Take the expt data and the baseline data and merge the two together
#'    warning about any missing genes
#'    
#' @param sample_info data.frame - sample info
#' @param expt_count_data data.frame - expt count data
#' @param baseline_data SummarizedExperiment object - baseline data
#' @param session session_object
#' 
#' @return SummarizedExperiment
#'
#' @examples
#' merge_with_baseline( sample_info, expt_count_data, baseline_data, session_obj )
#' 
#' @export
#'
merge_with_baseline <- function( sample_info, expt_count_data, baseline_data, session_obj ) {
  
  common_genes <- intersect(rownames(data), rownames(Mm_baseline))
  
}