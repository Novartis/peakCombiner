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
test_data_filtered <- peakCombiner::filter_regions(
  data = test_data_center_expand,
  include_by_chromosome_name = NULL,
  exclude_by_blacklist = "hg38",
  include_above_score_cutoff = NULL,
  include_top_n_scoring = NULL
)
##
result_colnames <- colnames(test_data_filtered)
##
test_data_combined <- peakCombiner::combine_regions(
  data = test_data_filtered,
  found_in_samples = 2,
  combined_center = "nearest"
)
##
test_data_combined_ce <- peakCombiner::center_expand_regions(
  data = test_data_combined,
  center_by = "center_column",
  expand_by = NULL
)
##
### -----------------------------------------------------------------------###
### Test input
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
  expect_equal(round(mean(data$center), 0), 2458)
  expect_identical(nrow(data), 52L)
  expect_identical(data$start[1], 350L)
  ##
  test_counts_left <- test_data_filtered |>
    dplyr::group_by(sample_name) |>
    dplyr::summarise(counts = dplyr::n()) |>
    dplyr::filter(sample_name == "treatment_rep1") |>
    dplyr::pull(counts)
  expect_identical(test_counts_left, 9L)
})
##
### --------------------------------------------------------------------------###
##
test_that("Output data frame is correct for data_prepared", {
  ##
  data <- test_data_prepared
  ##
  expect_no_error(peakCombiner::filter_regions(
    data = data,
    exclude_by_blacklist = "hg38",
    include_by_chromosome_name = NULL,
    include_above_score_cutoff = NULL,
    include_top_n_scoring = NULL
  ))
  ##
  result <- peakCombiner::filter_regions(
    data = data,
    exclude_by_blacklist = "hg38",
    include_by_chromosome_name = NULL,
    include_above_score_cutoff = NULL,
    include_top_n_scoring = NULL
  )
  ##
  expect_identical(nrow(result), 52L)
  expect_identical(result$start[9], 301L)
})
##
test_that("Output data frame is correct for data_center_expand", {
  ##
  data <- test_data_center_expand
  ##
  result <- peakCombiner::filter_regions(
    data = data,
    exclude_by_blacklist = "hg38",
    include_by_chromosome_name = NULL,
    include_above_score_cutoff = NULL,
    include_top_n_scoring = NULL
  )
  ##
  expect_identical(nrow(result), 52L)
  expect_identical(result$start[9], 250L)
})
##
test_that("Output data frame is correct for data_filtered", {
  ##
  data <- test_data_filtered
  result <- peakCombiner::filter_regions(
    data = data,
    exclude_by_blacklist = "hg38",
    include_by_chromosome_name = NULL,
    include_above_score_cutoff = NULL,
    include_top_n_scoring = NULL
  )
  ##
  expect_identical(nrow(result), 52L)
  expect_identical(result$start[9], 250L)
})
##
test_that("Output data frame is correct for data_combined", {
  ##
  data <- test_data_combined
  result <- peakCombiner::filter_regions(
    data = data,
    exclude_by_blacklist = "hg38",
    include_by_chromosome_name = NULL,
    include_above_score_cutoff = NULL,
    include_top_n_scoring = NULL
  )
  ##
  expect_identical(nrow(result), 10L)
  expect_identical(result$start[9], 250L)
})
##
test_that("Output data frame is correct for data_combined_ce", {
  ##
  data <- test_data_combined_ce
  result <- peakCombiner::filter_regions(
    data = data,
    exclude_by_blacklist = "hg38",
    include_by_chromosome_name = NULL,
    include_above_score_cutoff = NULL,
    include_top_n_scoring = NULL
  )
  ##
  expect_identical(nrow(result), 10L)
  expect_identical(result$start[9], 250L)
})
##
### -----------------------------------------------------------------------###
