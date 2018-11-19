context('Load Data')
library('rprojroot')
rootPath <- find_root(is_rstudio_project)

source(file.path(rootPath, 'R/load_data.R'))
sample_file <- file.path(rootPath, 'data/test-brd2-samples.txt')
count_file <- file.path(rootPath, 'data/test-brd2-counts.tsv')
session <- list(userData = list(testing = TRUE, debug = FALSE))

sample_info <- load_sample_data( sample_file, session )
test_that("load_sample_data", {
  expect_equal(class(sample_info), 'data.frame')
  expect_equal(dim(sample_info), c(8, 4))
  expect_equal(levels(sample_info$condition), c("hom", "wt"))
  expect_equal(levels(sample_info$sex), c("F", "M"))
  expect_equal(levels(sample_info$stage), paste(c("15", "20", "13", "19", "26", "22", "28"), "somites"))
})

count_data <- load_count_data(count_file, session)
test_that("load_count_data", {
  expect_equal(class(count_data), 'data.frame')
  expect_equal(dim(count_data), c(1999, 9))
})

expt_data <- load_data( sample_file, count_file, session )
test_that("load_data", {
  expect_equal(class(expt_data)[1], 'SummarizedExperiment')
  expect_equal(dim(expt_data), c(1999, 8))
  expect_equal(dim(colData(expt_data)), c(8, 4))
})

load(file.path(rootPath, 'data/test-baseline.rda'))
test_that("merge_with_baseline", {
  expect_warning(merge_with_baseline(expt_data, Mm_baseline_test, session),
                 "The following genes are present in the [A-Za-z]+ data, but not in the [A-Za-z]+ data:")
})

