##
### -----------------------------------------------------------------------###
### Prepare data for testing
### -----------------------------------------------------------------------###
## tweak the prepare_input_regions() function and re-load it
devtools::load_all()
##
### -----------------------------------------------------------------------###
### Prepare data for testing
### -----------------------------------------------------------------------###
##
test_data <- peakCombiner::syn_sample_sheet
samplesheet_colnames <- colnames(test_data)
##
all_colnames <- c(
  "chrom", "start", "end", "score", "strand", "summit", "sample_name"
)
##
data_prepared <- load_input_regions(
  data = test_data
)
##
### -----------------------------------------------------------------------###
### Test input
### -----------------------------------------------------------------------###
##
test_that("Test if function works with correct input", {
  expect_no_error(load_input_regions(
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
  expect_identical(names(test_data), input_colnames)
})
##
test_that("Input column 'sample_name' is a class 'character'.", {
  expect_true(is.character(test_data$sample_name))
})
##
test_that("Input column 'file_path' is a class 'character'", {
  expect_true(is.character(test_data$file_path))
})
##
test_that("Input column 'file_format' is a class 'character'", {
  expect_true(is.character(test_data$file_format))
})
##
### -----------------------------------------------------------------------###
##
test_that("Error occurs when 'data' does not exist.", {
  expect_error(load_input_regions(
    data = "nonexisting"
  ),)
})
##
test_that("Error occurs when 'data' has the wrong structure.", {
  expect_error(load_input_regions(
    data = tibble(1:10)
  ))
})
##
test_that("Error occurs when 'data' is a vector.", {
  expect_error(load_input_regions(
    data = as.vector(1:10)),)
})
##
test_that("Error occurs when 'data' is 'NULL'.", {
  expect_error(load_input_regions(
    data = NULL),)
})
##
test_that("Error occurs when 'data' is 'NA'.", {
  expect_error(load_input_regions(
    data = NA),)
})
##
### -----------------------------------------------------------------------###
### Test output
### -----------------------------------------------------------------------###
##
test_that("Column names of output data are identical with required once.", {
  expect_setequal(colnames(data_prepared), all_colnames)
})
##
test_that("Output data has the right number of columns", {
  expect_equal(ncol(data_prepared), 7)
  ##
})
##
test_that("Output data has the right class.", {
  expect_identical(class(data_prepared)[2], "tbl")
  ##
})
##
test_that("Output data has in column 'score', row 1 the correct value.", {
  expect_identical(round(data_prepared$score[1],0), 4)
})
##
test_that("Output data has the correct number of rows.", {
  expect_identical(nrow(data_prepared), 55L)
})
