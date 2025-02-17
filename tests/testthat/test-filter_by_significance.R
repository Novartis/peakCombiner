##
### -----------------------------------------------------------------------###
### Prepare data for testing
### -----------------------------------------------------------------------###
##
set.seed(1234)
##
required_colnames <- c(
  "chrom", "start", "end", "name", "score", "strand",
  "center", "sample_name"
)
##
data(syn_data_tibble, package = "peakCombiner")
test_data <- syn_data_tibble
input_colnames <- colnames(test_data)
##
test_data_prepared <- peakCombiner::prepare_input_regions(
  data = test_data
)
test_data_center_expand <- peakCombiner::center_expand_regions(
  data = test_data_prepared,
  center_by = "center_column",
  expand_by = NULL
)
##
input_colnames <- colnames(test_data_center_expand)
##
filter_by_significance <- 40
##
test_data_filtered <- peakCombiner:::filter_by_significance(
  data = test_data_center_expand,
  include_above_score_cutoff = filter_by_significance
)
##
result_colnames <- colnames(test_data_filtered)
##
### -----------------------------------------------------------------------###
### Test input
### -----------------------------------------------------------------------###
##
test_that("Test if function works with correct input", {
  expect_no_error(peakCombiner:::filter_by_significance(
    data = test_data_center_expand,
    include_above_score_cutoff = filter_by_significance
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
test_that("Required parameter 'filter_by_significance' has expected structure", {
  expect_no_error(peakCombiner:::filter_by_significance(
    data = test_data_filtered,
    include_above_score_cutoff = NULL
  ))
  expect_no_error(peakCombiner:::filter_by_significance(
    data = test_data_filtered,
    include_above_score_cutoff = 0
  ))
  ##
  expect_error(peakCombiner:::filter_by_significance(
    data = test_data_filtered,
    include_above_score_cutoff = NA
  ))
  expect_error(peakCombiner:::filter_by_significance(
    data = test_data_filtered,
    include_above_score_cutoff = "nonexisting"
  ))
  expect_error(peakCombiner:::filter_by_significance(
    data = test_data_filtered,
    include_above_score_cutoff = c(1, 2, 3)
  ))
  ##
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
  expect_equal(round(mean(data$center), 0), 2547)
  expect_identical(nrow(data), 38L)
  expect_identical(data$start[1], 4550)
})
##
### -----------------------------------------------------------------------###
