context('Helper Functions')
library('rprojroot')
rootPath <- find_root(is_rstudio_project)

source(file.path(rootPath, 'R/helper_functions.R'))

# create test character vectors
test_char_vec <- paste0(13:19, 'somites')
test_factor_vec_small <- factor(test_char_vec)
test_factor_vec_large <- factor(paste0(10:25, 'somites'))

return_value_small <- c( 'blue' = rgb(0,0.45,0.7),
                         'yellow' = rgb(0.95, 0.9, 0.25),
                         'vermillion' = rgb(0.8, 0.4, 0),
                         'purple' = rgb(0.8, 0.6, 0.7),
                         'blue_green' = rgb(0, 0.6, 0.5),
                         'sky_blue' = rgb(0.35, 0.7, 0.9),
                         'black' = rgb(0, 0, 0))
names(return_value_small) <- test_factor_vec_small
test_that("colour_palette", {
  expect_error(colour_palette(test_char_vec), 'Not a factor')
  expect_equal(length(colour_palette(test_factor_vec_small)), 7)
  expect_equal(colour_palette(test_factor_vec_small), return_value_small)
  expect_equal(length(colour_palette(test_factor_vec_large)), 16)
})

test_factor_shape <- factor(rep(c('wt', 'het', 'hom'), each = 3 ))
test_factor_shape_large <- factor(rep(c(letters[1:5], 'baseline'), each = 3 ))
test_that("shape_palette", {
  expect_error(shape_palette(test_char_vec), 'Not a factor')
  expect_error(shape_palette(test_factor_shape), 'None of the levels are "baseline"')
  expect_error(shape_palette(test_factor_shape_large), 'There are two many levels')
  test_factor_shape <- factor(c(as.character(test_factor_shape), rep('baseline', each = 3)))
  expect_equal(length(shape_palette(test_factor_shape)), 4)
  expect_equal(shape_palette(test_factor_shape), c('baseline' = 21, 'het' = 22, 'hom' = 23, 'wt' = 24 ))
  test_factor_shape <- factor(c(as.character(test_factor_shape), rep('baseline', each = 3)),
                              levels = c('wt', 'het', 'hom', 'baseline'))
  expect_equal(shape_palette(test_factor_shape), c('baseline' = 21, 'wt' = 22, 'het' = 23, 'hom' = 24))
})
