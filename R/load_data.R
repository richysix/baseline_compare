#' Load data from sample and count files
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
#' load_data( 'samples.txt', 'all.tsv', session_obj )
#'
#' @export
#'
load_data <- function( sample_file, count_file, session ){
  # read in sample info
  sample_info <- load_sample_data(sample_file, session)
  ## TO DO
  # make sure stage column is compatible with baseline data
  
  # load count data
  expt_count_data <- load_count_data(count_file, session)
  ## TO DO
  # check for missing genes and warn
  # make a SummarizedExperiment object
  counts <-
    expt_count_data[, grepl(' count$', colnames(expt_count_data)) &
                 !grepl('normalised count$', colnames(expt_count_data))]
  colnames(counts) <-
    gsub(' count$', '', colnames(counts))
  rownames(counts) <- expt_count_data[, 'Gene.ID']
  row_data <- expt_count_data[ , !grepl('count$', colnames(expt_count_data)), drop = FALSE ]
  
  expt_data <- SummarizedExperiment(
    assays = list(counts = as.matrix(counts)),
    rowData = DataFrame(row_data),
    colData = DataFrame(sample_info)
  )
  return(expt_data)
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
  }
  
  if (session$userData[['debug']]) {
    print('Function: load_sample_data')
    print(head(sample_info))
  }
  
  return(sample_info)
}

#' Check condition column
#'
#' \code{valid_condition_column} takes a vector and checks it is a valid condition column
#'
#'    The values are checked against the list of allowed values.
#'    
#' @param condition_column   factor/character - vector to check
#' 
#' @return logical TRUE/FALSE for whether vector passes
#'
#' @examples
#' valid_condition_column( sample_info$condition )
#'
#' @export
#'
valid_condition_column <- 
  function(condition_column, 
           allowed_values = c('hom', 'het', 'wt', 'mut', 'sib') ){
  return(all(condition_column %in% allowed_values))
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
#' @param expt_data SummarizedExperiment object - expt data
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
merge_with_baseline <- function( expt_data, baseline_data, session_obj ) {
  # find genes in common between expt data and baseline
  common_genes <- intersect(rownames(expt_data), rownames(baseline_data))
  expt_only_genes <- setdiff(rownames(expt_data), rownames(baseline_data))
  baseline_only <- setdiff(rownames(baseline_data), rownames(expt_data))
  
  # create warning for any missing genes
  msg <- NULL
  if ( length(expt_only_genes) > 0 ) {
    msg <- paste0('The following genes are present in the experimental data, but not in the Baseline data: ',
                  paste0(expt_only_genes, collapse = ', '))
    warning(msg)
  }
  if ( length(baseline_only) > 0 ) {
    msg <- paste0('The following genes are present in the Baseline data, but not in the experimental data: ',
                  paste0(baseline_only, collapse = ', '))
    warning(msg)
  }

  # subset each to common genes
  baseline_subset <- baseline_data[ common_genes, ]
  expt_subset <- expt_data[ common_genes, ]
  
  # combine count data
  counts <- cbind(assays(expt_subset)$counts, assays(baseline_subset)$counts)
  
  # colData - find common columns and rbind
  common_columns <- intersect(names(colData(expt_data)), names(colData(baseline_data)))
  col_data <- rbind( colData(expt_data)[ , common_columns ],
                     colData(baseline_data)[ , common_columns ] )
  # set levels of stage
  # col_data$stage <- factor(col_data$stage,
  #                          levels = levels( colData(baseline_data)$stage ) )

  # create new SummarizedExperiment object and return it
  merged_data <- 
    SummarizedExperiment(
      assays = list(counts = counts),
      rowData = rowData(baseline_subset),
      colData = col_data
    )
  
  return(merged_data)
}
