##
### -----------------------------------------------------------------------###
### Prepare data for testing
### -----------------------------------------------------------------------###
## tweak the prepare_input_regions() function and re-load it
#devtools::load_all()
library("tidyverse")
library("GenomicRanges")
set.seed(1234)
##
### -----------------------------------------------------------------------###
### Prepare data for testing
### -----------------------------------------------------------------------###

colnames_preloaded_df <- c(
  "chrom", "start", "end", "name", "score", "strand",
  "center", "sample_name"
)

allowed_file_format <- c("narrowpeak", "broadpeak", "bed")

data(syn_data_bed)
samplesheet_test <- syn_data_bed

test_sample_sheet <- prepare_input_regions(
  data = samplesheet_test[1, ]
)

data(syn_data_tibble)
test_data <- syn_data_tibble
input_colnames <- colnames(test_data)

test_data_prepared <- prepare_input_regions(
  data = test_data
)

restult_colnames <- colnames(test_data_prepared)

### -----------------------------------------------------------------------###
### Test input
### -----------------------------------------------------------------------###
### Test pre-loaded data frame
test_that("Test if function works with correct input", {
  expect_no_error(prepare_input_regions(
    data = test_data
  ))
})

test_that("Input data has at least 8 number of columns", {
  expect_equal(length(colnames(test_data)), 8)
})

test_that("Column names of input data are identical with required once.", {
  expect_true(all(colnames_preloaded_df %in% names(test_data)))
})

### -----------------------------------------------------------------------###
### Test pre-loaded gRanges
### -----------------------------------------------------------------------###

test_that("Input data has the right number of columns", {
  expect_equal(length(input_colnames), 8)
})

test_that("Input column 'chr' is a class 'character'.", {
  expect_true(is.character(test_data$chrom))
})

test_that("Input column 'start' is a class 'numeric'.", {
  expect_true(is.numeric(test_data$start))
})

test_that("Input column 'end' is a class 'numeric'.", {
  expect_true(is.numeric(test_data$end))
})

test_that("Input column 'score' is a class 'numeric'.", {
  expect_true(is.numeric(test_data$score))
})

test_that("Input column 'strand' is a class 'character'.", {
  expect_true(is.character(test_data$strand))
})

test_that("Input column 'center' is a class 'numeric'.", {
  expect_true(is.numeric(test_data$center))
})

test_that("Input column 'sample_name' is a class 'character'.", {
  expect_true(is.character(test_data$sample_name))
})

### -----------------------------------------------------------------------###
### Test output
### -----------------------------------------------------------------------###

test_that("Output data frame has the correct structure.", {
  expect_no_error(check_data_structure(test_data_prepared))
})

test_that("Column names of output data are identical with required once.", {
  expect_setequal(colnames(test_data_prepared), restult_colnames)
})

test_that("Output data has the right number of columns", {
  expect_equal(ncol(test_data_prepared), 8)
})

test_that("Output data has the correct class", {
  expect_identical(class(test_data_prepared)[2], "tbl")
})

test_that("Ouput column 'chrom' is a class 'character'.", {
  expect_true(is.character(test_data_prepared$chrom))
})

test_that("Ouput column 'start' is a class 'numeric'.", {
  expect_true(is.numeric(test_data_prepared$start))
})

test_that("Ouput column 'end' is a class 'numeric'.", {
  expect_true(is.numeric(test_data_prepared$end))
})

test_that("Ouput column 'name' is a class 'character'.", {
  expect_true(is.character(test_data_prepared$name))
})

test_that("Ouput column 'score' is a class 'numeric'.", {
  expect_true(is.numeric(test_data_prepared$score))
})

test_that("Ouput column 'strand' is a class 'character'.", {
  expect_true(is.character(test_data_prepared$strand))
})

test_that("Ouput column 'center' is a class 'numeric'.", {
  expect_true(is.numeric(test_data_prepared$center))
})

test_that("Ouput column 'sample_name' is a class 'character'.", {
  expect_true(is.character(test_data_prepared$sample_name))
})

test_that("The mean of all output centers.", {
  expect_equal(round(mean(test_data_prepared$center), 0), 2452)
})

test_that("The number of rows in the output file.", {
  expect_identical(nrow(test_data_prepared), 52L)
})

### -----------------------------------------------------------------------###
