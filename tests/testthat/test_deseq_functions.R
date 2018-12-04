context('DESeq Functions')
library('rprojroot')
rootPath <- find_root(is_rstudio_project)

source(file.path(rootPath, 'R/deseq_functions.R'))

session_obj <- list(userData = list(testing = TRUE, debug = FALSE))
load(file.path(rootPath, 'data/test-brd2-data.rda'))
expt_only_dds <- create_new_DESeq2DataSet(expt_data, baseline_data = NULL,
                                          gender_column = NULL, groups = NULL,
                                          condition_column = 'condition',
                                          session_obj = session_obj)
test_that("create from expt data only", {
  expect_s4_class(expt_only_dds, class = 'DESeqDataSet')
  expect_equal(dim(colData(expt_only_dds)), c(8, 4), info = 'colData')
  expect_equal(dim(rowData(expt_only_dds)), c(1999, 1), info = 'rowData')
  expect_equal(dim(counts(expt_only_dds)), c(1999, 8), info = 'counts')
  expect_equal(design(expt_only_dds), as.formula('~ condition'), info = 'design formula')
})

# rerun with gender_column specified to match other datasets
expt_only_dds <- create_new_DESeq2DataSet(expt_data, baseline_data = NULL,
                                          gender_column = 'sex', groups = NULL,
                                          condition_column = 'condition',
                                          session_obj = session_obj)


load(file.path(rootPath, 'data/Mm_baseline_data.rda'))
# subset expt_data and Mm_baseline to common_genes to avoid warning
common_genes <- intersect(rownames(expt_data), rownames(Mm_baseline))
baseline_subset <- Mm_baseline[ common_genes, ]
expt_subset <- expt_data[ common_genes, ]

expt_plus_baseline_dds <- 
  create_new_DESeq2DataSet(expt_subset, baseline_data = baseline_subset,
                          gender_column = 'sex', groups = c('sex'),
                          condition_column = 'condition', session_obj = session_obj)
test_that("create from expt_data plus baseline", {
  expect_s4_class(object = expt_plus_baseline_dds, class = 'DESeqDataSet')
  expect_equal(dim(colData(expt_plus_baseline_dds)), c(36, 4), info = 'colData')
  expect_equal(dim(rowData(expt_plus_baseline_dds)), c(1999, 8), info = 'rowData')
  expect_equal(dim(counts(expt_plus_baseline_dds)), c(1999, 36), info = 'counts')
  expect_equal(design(expt_plus_baseline_dds), 
               as.formula('~ sex + condition'), info = 'design formula')
  expect_warning(create_new_DESeq2DataSet(expt_subset, baseline_data = Mm_baseline,
                                          gender_column = 'sex', groups = c('sex'),
                                          condition_column = 'condition', 
                                          session_obj = session_obj),
                 "The following genes are present in the [A-Za-z]+ data, but not in the [A-Za-z]+ data:")
})

expt_plus_baseline_with_stage_dds <- 
  create_new_DESeq2DataSet(expt_subset, baseline_data = baseline_subset,
                           gender_column = 'sex', groups = c('sex', 'stage'),
                           condition_column = 'condition', session_obj = session_obj)
test_that("create from expt_data plus baseline with stage", {
  expect_s4_class(object = expt_plus_baseline_with_stage_dds, class = 'DESeqDataSet')
  expect_equal(dim(colData(expt_plus_baseline_with_stage_dds)), c(36, 4), info = 'colData')
  expect_equal(dim(rowData(expt_plus_baseline_with_stage_dds)), c(1999, 8), info = 'rowData')
  expect_equal(dim(counts(expt_plus_baseline_with_stage_dds)), c(1999, 36), info = 'counts')
  expect_equal(design(expt_plus_baseline_with_stage_dds), 
               as.formula('~ sex + stage + condition'), info = 'design formula')
  expect_warning(create_new_DESeq2DataSet(expt_subset, baseline_data = Mm_baseline,
                                          gender_column = 'sex', groups = c('sex', 'stage'),
                                          condition_column = 'condition', session_obj = session_obj),
                 "The following genes are present in the [A-Za-z]+ data, but not in the [A-Za-z]+ data:")
})

expt_plus_all_baseline_dds <- 
  create_new_DESeq2DataSet(expt_subset, baseline_data = baseline_subset,
                           gender_column = NULL, groups = NULL,
                           condition_column = 'condition', 
                           match_stages = FALSE, session_obj = session_obj)
expt_plus_baseline_matched_stage_dds <- 
  create_new_DESeq2DataSet(expt_subset, baseline_data = baseline_subset,
                         gender_column = NULL, groups = NULL,
                         condition_column = 'condition', 
                         match_stages = TRUE, session_obj = session_obj)
test_that("create from expt_data plus baseline, match or not stages", {
  expect_s4_class(object = expt_plus_all_baseline_dds, class = 'DESeqDataSet')
  expect_equal(dim(colData(expt_plus_all_baseline_dds)), c(119, 4), info = 'colData')
  expect_equal(levels(colData(expt_plus_all_baseline_dds)[['stage']]), 
               levels(colData(Mm_baseline)[['stage']]), info = 'colData stage levels')
  expect_equal(dim(rowData(expt_plus_all_baseline_dds)), c(1999, 8), info = 'rowData')
  expect_equal(dim(counts(expt_plus_all_baseline_dds)), c(1999, 119), info = 'counts')
  expect_equal(dim(colData(expt_plus_baseline_matched_stage_dds)),
               c(36, 4), info = 'match stages')
  expect_equal(levels(colData(expt_plus_baseline_matched_stage_dds)[['stage']]), 
                paste0(c(13, 15, 19, 20, 22, 26, 28), 'somites'), info = 'match stages levels')
})

expt_data_missing_stage <- expt_data
colData(expt_data_missing_stage)[['stage']] <- 
  factor(sub("15", "1", colData(expt_data_missing_stage)[['stage']]),
        levels = sub("15", "1", levels(colData(expt_data_missing_stage)[['stage']]) ) )

test_that("missing stage in expt data", {
  expect_warning(
    create_new_DESeq2DataSet(expt_data_missing_stage, baseline_data = baseline_subset,
                              gender_column = NULL, groups = NULL,
                              condition_column = 'condition', 
                              match_stages = TRUE, session_obj = session_obj),
    'The following stages are present in the [A-Za-z]+ data, but not in the [A-Za-z]+ data:', 
    info = 'missing stage'
  )
})

test_that("Model matrix not full rank", {
  expt_data_not_full_rank <- expt_data
  colData(expt_data_not_full_rank)[['stage']] <- 
    factor(c(rep('15somites', 5), rep('25somites', 3)), 
           levels = c('15somites', '25somites'))
  expect_warning(create_new_DESeq2DataSet(expt_data_not_full_rank, 
                                          gender_column = 'sex', 
                                          groups = c('sex', 'stage'),
                                          condition_column = 'condition',
                                          session_obj = session_obj),
                 "Model matrix not full rank"
  )
})

# Test run_deseq
deseq_res <- run_deseq(expt_only_dds, 'mut', 'sib', 0.05)
test_that("DESeq results", {
  expect_equal(sum(deseq_res$result$padj < 0.05 & !is.na(deseq_res$result$padj)),
               524)
  expect_equal(head( order(deseq_res$result$padj) ),
               c(130, 1938, 733, 1294, 1384, 599) )
})

# test overlap
deseq_datasets <- list(
  expt_only_dds = expt_only_dds,
  expt_plus_baseline_dds = expt_plus_baseline_dds,
  expt_plus_baseline_with_stage_dds = expt_plus_baseline_with_stage_dds
)
overlapped_results <- overlap_deseq_results(deseq_datasets, 'mut', 'sib', 0.05, session_obj)

test_that("DESeq2 overlaps", {
  expect_equal(class(overlapped_results), 'list')
  expect_equal(length(overlapped_results), 5)
  expect_equal(names(overlapped_results), 
               c('expt_only_res', 'plus_baseline_res', 'plus_baseline_with_stage_res',
                 'overlaps', 'results_tables'))
  # check overlaps
  expect_equal(names(overlapped_results[['overlaps']]),
                     c('mutant_response', 'delay', 'no_delay', 'discard'))
  expect_equal(length(overlapped_results[['overlaps']][['mutant_response']]), 162)
  expect_equal(length(overlapped_results[['overlaps']][['delay']]), 74)
  expect_equal(length(overlapped_results[['overlaps']][['no_delay']]), 172)
  expect_equal(length(overlapped_results[['overlaps']][['discard']]), 235)
  # check results tables
  expect_equal(length(overlapped_results[['results_tables']]), 6)
  expect_equal(dim(overlapped_results[['results_tables']][['unprocessed']]), c(524, 8))
  expect_equal(dim(overlapped_results[['results_tables']][['mutant_response']]), c(162, 12))
  expect_equal(dim(overlapped_results[['results_tables']][['delay']]), c(74, 12))
  expect_equal(dim(overlapped_results[['results_tables']][['no_delay']]), c(172, 12))
  expect_equal(dim(overlapped_results[['results_tables']][['discard']]), c(235, 12))
  expect_equal(dim(overlapped_results[['results_tables']][['all_genes']]), c(1999, 12))
})
