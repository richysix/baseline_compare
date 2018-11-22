library('DESeq2')

#' Create a DESeq2 DataSet from the supplied data
#'
#' \code{create_new_DESeq2DataSet}
#'
#'    This function takes the experimental and combined datasets and runs DESeq2 three
#'    different ways on the data. 
#'    1. The experimental samples on their own  (model: either ~ condition or ~ sex + condition)
#'    2. The experimental and baseline samples  (model: either ~ condition or ~ sex + condition)
#'    3. The experimental and baseline samples  (model: either ~ stage + condition or ~ sex + stage + condition)
#'    
#' @param expt_data        SummarizedExperiment - object containing the data for experiment samples only
#' @param baseline_data    SummarizedExperiment - object containing the data for the baseline data
#' @param gender_column    character - name of the column representing the gender of the experimental samples if present
#'                         This is NULL if gender should not be used.
#' @param groups           vector - variables to use in the design formula
#' @param condition_column character - name of the columns in the experimental sample information to use as condition
#' @param match_stage      logical - wether to only include baseline samples that match the stages in the experimental data
#' @param session_obj      session_object
#' 
#' @return A DESeq2DataSet object
#'
#' @examples
#' create_new_DESeq2DataSet( expt_data, baseline_data = NULL, gender_column = NULL, groups = NULL, condition_column = 'condition', session_obj )
#' 
#' create_new_DESeq2DataSet( expt_data, baseline_data = Mm_baseline, gender_column = 'sex', groups = NULL, condition_column = 'treatment', session_obj )
#' 
#' create_new_DESeq2DataSet( expt_data, baseline_data = Mm_baseline, gender_column = 'gender', groups = c('gender', 'stage'), condition_column = 'condition', session_obj )
#'
create_new_DESeq2DataSet <- function( expt_data, baseline_data = NULL,
                                      gender_column = NULL, groups = NULL,
                                      condition_column = 'condition', 
                                      match_stages = TRUE, session_obj ) {
  # change names of expt data sample info to standardize
  if (!is.null(gender_column)) {
    names(colData(expt_data))[ names(colData(expt_data)) == gender_column ] <- 'sex'
    if (is.null(groups)) {
      groups <- c('sex')
    }
  }
  if (!is.null(condition_column)) {
    names(colData(expt_data))[ names(colData(expt_data)) == condition_column ] <- 'condition'
  }
  
  # create design formula
  if (!is.null(groups)) {
    design_formula <- as.formula(paste0('~', paste0(groups, collapse = ' + '), ' + condition') )
  } else {
    design_formula <- as.formula('~ condition')
  }
  
  if (!is.null(baseline_data)) {
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
    
    # check for stages missing from the baseline data
    expt_only_stages <- setdiff(levels(colData(expt_data)[['stage']]), 
                                levels(colData(baseline_data)[['stage']]) )
    if ( length(expt_only_stages) > 0 ) {
      msg <- paste0('The following stages are present in the experimental data, but not in the Baseline data: ',
                    paste0(expt_only_stages, collapse = ', '))
      warning(msg)
    }
    
    # find baseline samples that match the stages in expt_data
    if (match_stages) {
      samples_to_keep <- 
        Reduce('|', lapply(levels(colData(expt_data)[['stage']]), 
                    function(stage, baseline_subset){ 
                      colData(baseline_subset)[['stage']] == stage 
                    }, 
              baseline_subset
             ))
      baseline_subset <- baseline_subset[ , samples_to_keep ]
      colData(baseline_subset) <- droplevels(colData(baseline_subset))
    }
    
    # combine count data
    counts <- cbind(assays(expt_subset)$counts, assays(baseline_subset)$counts)
    
    # colData - find common columns and rbind
    common_columns <- intersect(names(colData(expt_data)), names(colData(baseline_subset)))
    col_data <- rbind( colData(expt_data)[ , common_columns ],
                       colData(baseline_subset)[ , common_columns ] )
    # set levels to match the order of baseline
    if ( length(expt_only_stages) == 0 ) {
      col_data$stage <- factor(col_data$stage, levels = levels(colData(baseline_subset)[['stage']]))
    } else {
      all_levels <- unique(col_data$stage)
      stage_numbers <- as.numeric(sub("somites", "", all_levels))
      col_data$stage <- factor(col_data$stage, 
                               levels = all_levels[ order(stage_numbers) ])
    }
    
    row_data <- rowData(baseline_subset)
    se <- SummarizedExperiment(
      assays = list(counts = counts),
      rowData = row_data,
      colData = col_data
    )
  } else {
    se <- expt_data
  }
  
  dds <- tryCatch( DESeqDataSet(se, design = design_formula),
                   error = function(err){
                     if( grepl('the model matrix is not full rank', err$message ) ){
                       if (design_formula == '~condition') {
                         stop(err$message)
                       }
                       warning("Model matrix not full rank. Running with just condition...\n")
                       design_formula <- as.formula('~condition')
                       DESeqDataSet(se, design = design_formula)
                     } else {
                       stop(err$message)
                     }
                   }
  )
                 
  return(dds)
}

#' Run DESeq2
#'
#' \code{run_deseq} Runs DESeq2 on the supplied DESeqDataSet and returns the results
#'
#' @param dds DESeqDataSet object
#' @param expt_condition character - A level of the condition variable to compare
#' @param ctrl_condition character - A level of the condition variable to compare to
#' 
#' @return list containing:
#'         deseq = DESeqDataSet object
#'         results = DESeqResults object
#'
#' @examples
#' run_deseq( deseq_dataset, expt_condition, ctrl_condition )
#'
run_deseq <- function(dds, expt_condition, ctrl_condition) {
  dds <- tryCatch( DESeq(dds),
                   error = function(err){
                     if( grepl('the model matrix is not full rank', err$message ) ){
                       cat("Model matrix not full rank. Running with just condition...\n")
                       design(dds) <- formula('~ condition')
                       DESeq(dds)
                     } else {
                       stop(err$message)
                     }
                   }
  )
  
  # Write out results for specified pair of conditions
  res <- results(dds, contrast=c("condition", expt_condition, ctrl_condition))
  return(list(deseq = dds, result = res))
}
