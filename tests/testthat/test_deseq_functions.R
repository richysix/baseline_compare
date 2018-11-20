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

load(file.path(rootPath, 'data/test-baseline.rda'))
# subset expt_data and baseline_test to common_genes to avoid warning
common_genes <- intersect(rownames(expt_data), rownames(Mm_baseline_test))
baseline_subset <- Mm_baseline_test[ common_genes, ]
expt_subset <- expt_data[ common_genes, ]

expt_plus_baseline_dds <- 
  create_new_DESeq2DataSet(expt_subset, baseline_data = baseline_subset,
                          gender_column = 'sex', groups = c('sex'),
                          condition_column = 'condition', session_obj = session_obj)
test_that("create from expt_data plus baseline", {
  expect_s4_class(object = expt_plus_baseline_dds, class = 'DESeqDataSet')
  expect_equal(dim(colData(expt_plus_baseline_dds)), c(36, 4), info = 'colData')
  expect_equal(dim(rowData(expt_plus_baseline_dds)), c(1998, 8), info = 'rowData')
  expect_equal(dim(counts(expt_plus_baseline_dds)), c(1998, 36), info = 'counts')
  expect_equal(design(expt_plus_baseline_dds), 
               as.formula('~ sex + condition'), info = 'design formula')
  expect_warning(create_new_DESeq2DataSet(expt_subset, baseline_data = Mm_baseline_test,
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
  expect_equal(dim(rowData(expt_plus_baseline_with_stage_dds)), c(1998, 8), info = 'rowData')
  expect_equal(dim(counts(expt_plus_baseline_with_stage_dds)), c(1998, 36), info = 'counts')
  expect_equal(design(expt_plus_baseline_with_stage_dds), 
               as.formula('~ sex + stage + condition'), info = 'design formula')
  expect_warning(create_new_DESeq2DataSet(expt_subset, baseline_data = Mm_baseline_test,
                                          gender_column = 'sex', groups = c('sex', 'stage'),
                                          condition_column = 'condition', session_obj = session_obj),
                 "The following genes are present in the [A-Za-z]+ data, but not in the [A-Za-z]+ data:")
})

load(file.path(rootPath, 'data/Mm_baseline_data.rda'))
# subset Mm_baseline to common_genes to avoid warning
common_genes <- intersect(rownames(expt_subset), rownames(Mm_baseline))
Mm_baseline_subset <- Mm_baseline[ common_genes, ]

expt_plus_all_baseline_dds <- 
  create_new_DESeq2DataSet(expt_subset, baseline_data = Mm_baseline_subset,
                           gender_column = NULL, groups = NULL,
                           condition_column = 'condition', 
                           match_stages = FALSE, session_obj = session_obj)
expt_plus_baseline_matched_stage_dds <- 
  create_new_DESeq2DataSet(expt_subset, baseline_data = Mm_baseline_subset,
                         gender_column = NULL, groups = NULL,
                         condition_column = 'condition', 
                         match_stages = TRUE, session_obj = session_obj)
test_that("create from expt_data plus baseline, match or not stages", {
  expect_s4_class(object = expt_plus_all_baseline_dds, class = 'DESeqDataSet')
  expect_equal(dim(colData(expt_plus_all_baseline_dds)), c(119, 4), info = 'colData')
  expect_equal(dim(rowData(expt_plus_all_baseline_dds)), c(1998, 8), info = 'rowData')
  expect_equal(dim(counts(expt_plus_all_baseline_dds)), c(1998, 119), info = 'counts')
  expect_equal(dim(colData(expt_plus_baseline_matched_stage_dds)),
               c(36, 4), info = 'match stages')
})

