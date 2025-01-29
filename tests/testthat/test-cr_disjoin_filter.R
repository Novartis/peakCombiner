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
  "chrom", "start", "end", "name", "score", "strand", "center", "sample_name"
)
##
data(syn_data_tibble)
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
test_data_filtered <- peakCombiner::filter_regions(
  data = test_data_center_expand,
  exclude_by_blacklist = "hg38", # "hg38",
  include_by_chromosome_name = NULL,
  include_above_score_cutoff = NULL,
  include_top_n_scoring = NULL
)

##
test_data_disjoin_filter <- peakCombiner:::cr_disjoin_filter(
  data = test_data_filtered,
  found_in_samples = 2
)
##
result_colnames <- colnames(test_data_disjoin_filter)
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
test_that("Parameter 'found_in_samples' has the correct structure", {
  expect_no_error(peakCombiner:::cr_disjoin_filter(
    data = test_data_filtered,
    found_in_samples = 3
  ))
  expect_error(peakCombiner:::cr_disjoin_filter(
    data = test_data_filtered,
    found_in_samples = 0
  ), "Arg")
  expect_error(peakCombiner:::cr_disjoin_filter(
    data = test_data_filtered,
    found_in_samples = NULL
  ), "Arg")
  expect_error(peakCombiner:::cr_disjoin_filter(
    data = test_data_filtered,
    found_in_samples = NA
  ), )
  expect_error(peakCombiner:::cr_disjoin_filter(
    data = test_data_filtered,
    found_in_samples = c(1, 2, 3)
  ), "'")
  expect_error(peakCombiner:::cr_disjoin_filter(
    data = test_data_filtered,
    found_in_samples = test_data_filtered
  ), "Arg")
})
### -----------------------------------------------------------------------###
### Test Output
### -----------------------------------------------------------------------###
##
test_that("Output data frame is correct", {
  data <- test_data_disjoin_filter |>
    dplyr::mutate(chrom = as.character(chrom))
  ##
  expect_setequal(colnames(data), result_colnames)
  expect_equal(ncol(data), 12)
  ##
  expect_identical(class(data)[2], "tbl")
  ##
  expect_true(is.character(data$chrom))
  expect_true(is.numeric(data$start))
  expect_true(is.numeric(data$end))
  expect_true(is.character(data$sample_name))
  ##
  expect_identical(nrow(data), as.integer(113))
  expect_identical(data$start[1], 150)
  ##
  test_counts_left <- test_data_filtered |>
    dplyr::group_by(sample_name) |>
    dplyr::summarise(counts = dplyr::n()) |>
    dplyr::filter(sample_name == "treatment_rep1") |>
    dplyr::pull(counts)
  expect_identical(test_counts_left, as.integer(9))
})
##
### -----------------------------------------------------------------------###
