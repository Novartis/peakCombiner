##
### -----------------------------------------------------------------------###
### Prepare data for testing
### -----------------------------------------------------------------------###
##
set.seed(1234)
##
data(syn_data_tibble)
test_data <- syn_data_tibble
##
all_colnames <- c(
  "chrom", "start", "end", "name","score", "strand", "center", "sample_name"
)
input_colnames <- c(
  "chrom", "start", "end", "sample_name"
)
##
data_prepared <- test_data
##
### -----------------------------------------------------------------------###
### Test input
### -----------------------------------------------------------------------###
##
test_that("Test if function works with correct input", {
  expect_error(peakCombiner::load_input_regions(
    data = test_data
  ))
})
##
test_that("Parameter `score_colname` has a length of 1.", {
  expect_equal(length("score_colname"), 1)
})
##
test_that("Parameter `score_colname` is a character vector.", {
  expect_equal(class("score_colname"), "character")
})
##
test_that("Parameter `score_colname` is not numeric.", {
  expect_false(is.numeric("score_colname"))
})
##
### -----------------------------------------------------------------------###
##
test_that("Input data has exact three columns.", {
  expect_equal(length(input_colnames), 4)
})
##
test_that("Input data colnames are the expected once.", {
  expect_identical(names(test_data), all_colnames)
})
##
test_that("Input column 'sample_name' is a class 'character'.", {
  expect_true(is.character(test_data$sample_name))
})
##
### -----------------------------------------------------------------------###
##
test_that("Error occurs when 'data' does not exist.", {
  expect_error(peakCombiner::load_input_regions(
    data = "nonexisting"
  ), )
})
##
test_that("Error occurs when 'data' has the wrong structure.", {
  expect_error(peakCombiner::load_input_regions(
    data = tibble(1:10)
  ))
})
##
test_that("Error occurs when 'data' is a vector.", {
  expect_error(peakCombiner::load_input_regions(
    data = as.vector(1:10)
  ), )
})
##
test_that("Error occurs when 'data' is 'NULL'.", {
  expect_error(peakCombiner::load_input_regions(
    data = NULL
  ), )
})
##
test_that("Error occurs when 'data' is 'NA'.", {
  expect_error(peakCombiner::load_input_regions(
    data = NA
  ), )
})
##
### -----------------------------------------------------------------------###
### Test output
### -----------------------------------------------------------------------###
##
test_that("Output data has the right class.", {
  expect_identical(class(data_prepared)[2], "tbl")
  ##
})
##
test_that("Output data has in column 'score', row 1 the correct value.", {
  expect_identical(round(data_prepared$score[1], 0), 100)
})
##
test_that("Output data has the correct number of rows.", {
  expect_identical(nrow(data_prepared), 55L)
})

