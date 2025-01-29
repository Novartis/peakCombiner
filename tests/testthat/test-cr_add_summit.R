##
### -----------------------------------------------------------------------###
### Prepare data for testing
### -----------------------------------------------------------------------###
##
set.seed(1234)
##
input_colnames <- c(
  "chrom", "start", "end", "width", "strand", "name", "center", "score"
)
##
required_colnames <- c(
  "chrom", "start", "end", "name", "score", "strand",
  "center", "sample_name"
)
output_colnames <- c(
  "chrom", "start", "end", "strand", "name", "score", "center",
  "sample_name", "input_names"
)
##
data(syn_data_tibble, package = "peakCombiner")
test_data <- syn_data_tibble
##
test_data_prepared <- peakCombiner::prepare_input_regions(
  data = test_data
)

test_data_center_expand <- peakCombiner::center_expand_regions(
  data = test_data_prepared,
  center_by = "center_column",
  expand_by = NULL
)

test_data_filtered <- peakCombiner::filter_regions(
  data = test_data_center_expand,
  exclude_by_blacklist = "hg38",
  include_by_chromosome_name = c("chr1", "chr10", "chr2", "chr42"),
  include_above_score_cutoff = NULL,
  include_top_n_scoring = NULL
)

test_data_disjoin_filter <- peakCombiner:::cr_disjoin_filter(data = test_data_filtered, found_in_samples = 2)
test_data_reduce <- peakCombiner:::cr_reduce(data = test_data_disjoin_filter)
test_data_overlap <- peakCombiner:::cr_overlap_with_summits(
  data = test_data_reduce,
  input = test_data_filtered
)
##
test_data_combined_with_summit <- peakCombiner:::cr_add_summit(
  data = test_data_overlap,
  input = test_data_filtered,
  combined_center = "nearest",
  annotate_with_input_names = TRUE,
  combined_sample_name = "combined"
)
##
### -----------------------------------------------------------------------###
### Test input
### -----------------------------------------------------------------------###
##
test_that("Input data frame has the expected structure", {
  ##
  data <- test_data_overlap
  ##
  expect_equal(length(colnames(data)), 8)
  expect_identical(names(data), input_colnames)
  expect_true(is.factor(data$chrom))
  expect_true(is.numeric(data$start))
  expect_true(is.numeric(data$end))
  expect_true(is.character(data$name))
  expect_true(sum(stringr::str_detect(data$name, "|")) > 0)
  ##
})
##
test_that("Meta data frame has the expected structure", {
  ##
  data <- test_data_filtered
  ##
  expect_equal(length(colnames(data)), 8)
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
test_that("Parameter 'center' has the expected structure", {
  expect_no_error(peakCombiner:::cr_add_summit(
    data = test_data_overlap,
    input = test_data_filtered,
    combined_center = "STRONGEST"
  ))
  expect_no_error(peakCombiner:::cr_add_summit(
    data = test_data_overlap,
    input = test_data_filtered,
    combined_center = "middle"
  ))
  ##
  expect_error(peakCombiner:::cr_add_summit(
    data = test_data_overlap,
    input = test_data_filtered,
    combined_center = mean
  ))
  expect_error(peakCombiner:::cr_add_summit(
    data = test_data_overlap,
    input = test_data_filtered,
    combined_center = 2
  ))
  expect_error(peakCombiner:::cr_add_summit(
    data = test_data_overlap,
    input = test_data_filtered,
    combined_center = c(1, 2, 3)
  ), "`")
  expect_error(peakCombiner:::cr_add_summit(
    data = test_data_overlap,
    input = test_data_filtered,
    combined_center = NULL
  ), "`")
  expect_error(peakCombiner:::cr_add_summit(
    data = test_data_overlap,
    input = test_data_filtered,
    combined_center = NA
  ), "`")
})
##
### -----------------------------------------------------------------------###
### Test Output
### -----------------------------------------------------------------------###
##
test_that("Output data frame is correct", {
  ##
  data <- test_data_combined_with_summit
  ##
  expect_setequal(colnames(data), output_colnames)
  expect_equal(ncol(data), 9)
  ##
  expect_identical(class(data)[2], "tbl")
  ##
  expect_true(is.character(data$chrom))
  expect_true(is.numeric(data$start))
  expect_true(is.numeric(data$end))
  expect_true(is.character(data$input_names))
  ##
  expect_identical(nrow(data), as.integer(8))
  expect_identical(data$center[1], 500)
  expect_identical(round(sum(data$score), 0), 660)
  ##
})
##
test_that("Output data results with different summits", {
  data <- peakCombiner:::cr_add_summit(
    data = test_data_overlap,
    input = test_data_filtered,
    combined_center = "nearest"
  )
  expect_identical(data$center[7], 500)
  ##
  data <- peakCombiner:::cr_add_summit(
    data = test_data_overlap,
    input = test_data_filtered,
    combined_center = "strongest"
  )
  expect_identical(data$center[7], 600)
  ##
  data <- peakCombiner:::cr_add_summit(
    data = test_data_overlap,
    input = test_data_filtered,
    combined_center = "middle"
  )
  expect_identical(data$center[7], 550)
  ##
})
##
