##
### -----------------------------------------------------------------------###
### Prepare data for testing
### -----------------------------------------------------------------------###
##
library(peakCombiner)
##
set.seed(1234)
##
required_colnames <- c(
  "chrom", "start", "end", "name", "score", "strand",
  "center", "sample_name"
)
##
data(syn_data_tibble)
test_data <- syn_data_tibble
input_colnames <- colnames(test_data)
##
test_data_prepared <- peakCombiner::prepare_input_regions(
  data = test_data
)
##
test_data_center_expand <- peakCombiner::center_expand_regions(
  data = test_data_prepared,
  center_by = "center_column",
  expand_by = NULL
)
##
input_colnames <- colnames(test_data_center_expand)
##
keep_chromosomes <- c("chr1", "chr10", "chr42")
##
test_data_filtered <- peakCombiner:::filter_by_chromosome_names(
  data = test_data_center_expand,
  include_by_chromosome_name = keep_chromosomes
)
##
result_colnames <- colnames(test_data_filtered)
##
### -----------------------------------------------------------------------###
### Test input
### -----------------------------------------------------------------------###
##
test_that("Test if function works with correct input", {
  expect_no_error(peakCombiner:::filter_by_chromosome_names(
    data = test_data_center_expand,
    include_by_chromosome_name = keep_chromosomes
  ))
})
##
### -----------------------------------------------------------------------###
##
test_that("Input data frame has the expected structure", {
  data <- test_data_filtered
  ##
  expect_equal(length(input_colnames), 8)
  expect_identical(names(data), required_colnames)
  expect_true(is.character(data$chrom))
  expect_true(is.numeric(data$start))
  expect_true(is.numeric(data$end))
  expect_true(is.character(data$name))
  expect_true(is.numeric(data$score))
  expect_true(is.character(data$strand))
  expect_true(is.numeric(data$center))
  expect_true(is.character(data$sample_name))
  expect_true(sum(stringr::str_detect(data$name, "|")) > 0)
})
##
### -----------------------------------------------------------------------###
##
test_that("Required parameter 'filter_by_chromosome_names' has expected
          structure", {
  expect_no_error(peakCombiner:::filter_by_chromosome_names(
    data = test_data_filtered,
    include_by_chromosome_name = NULL
  ))
  expect_no_error(peakCombiner:::filter_by_chromosome_names(
    data = test_data_filtered,
    include_by_chromosome_name = "chr1"
  ))
  expect_no_error(peakCombiner:::filter_by_chromosome_names(
    data = test_data_filtered,
    include_by_chromosome_name = keep_chromosomes
  ))
  ##
  expect_error(peakCombiner:::filter_by_chromosome_names(
    data = test_data_filtered,
    include_by_chromosome_name = NA
  ))
})
##
### -----------------------------------------------------------------------###
### Test Output
### -----------------------------------------------------------------------###
##
test_that("Output data frame is correct", {
  data <- test_data_filtered
  ##
  expect_setequal(colnames(data), required_colnames)
  expect_equal(ncol(data), 8)
  ##
  expect_identical(class(data)[2], "tbl")
  ##
  expect_true(is.character(data$chrom))
  expect_true(is.numeric(data$start))
  expect_true(is.numeric(data$end))
  expect_true(is.character(data$name))
  expect_true(is.numeric(data$score))
  expect_true(is.character(data$strand))
  expect_true(is.numeric(data$center))
  expect_true(is.character(data$sample_name))
  ##
  expect_equal(round(mean(data$center), 0), 3168)
  expect_identical(nrow(data), 38L)
  expect_identical(data$start[1], 250)
})
##
### -----------------------------------------------------------------------###
