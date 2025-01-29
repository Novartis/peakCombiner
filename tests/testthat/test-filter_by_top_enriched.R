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
##
input_colnames <- colnames(test_data_center_expand)
##
test_data_filtered <- peakCombiner:::filter_by_top_enriched(
  data = test_data_center_expand,
  include_top_n_scoring = 10
)
##
result_colnames <- colnames(test_data_filtered)
##
table(test_data_center_expand$sample_name)
table(test_data_filtered$sample_name)
### -----------------------------------------------------------------------###
### Test input
### -----------------------------------------------------------------------###
##
test_that("Test if function works with correct input", {
  expect_no_error(peakCombiner:::filter_by_top_enriched(
    data = test_data_center_expand,
    include_top_n_scoring = 10
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
test_that("Required parameter 'filter_by_top_enriched' has expected structure", {
  expect_no_error(peakCombiner:::filter_by_top_enriched(
    data = test_data_center_expand,
    include_top_n_scoring = NULL
  ))
  expect_no_error(peakCombiner:::filter_by_top_enriched(
    data = test_data_center_expand,
    include_top_n_scoring = 5
  ))
  ##
  expect_error(peakCombiner:::filter_by_top_enriched(
    data = test_data_center_expand,
    include_top_n_scoring = 0
  ))
  expect_error(peakCombiner:::filter_by_top_enriched(
    data = test_data_center_expand,
    include_top_n_scoring = NA
  ))
  expect_error(peakCombiner:::filter_by_top_enriched(
    data = test_data_center_expand,
    include_top_n_scoring = "notexisting"
  ))
  expect_error(peakCombiner:::filter_by_top_enriched(
    data = test_data_center_expand,
    include_top_n_scoring = c(1, 2, 3)
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
  expect_equal(round(mean(data$center), 0), 2458)
  expect_identical(nrow(data), 52L)
  expect_identical(data$start[1], 350)
  ##
  test_counts_left <- test_data_filtered |>
    dplyr::group_by(sample_name) |>
    dplyr::summarise(counts = dplyr::n()) |>
    dplyr::filter(sample_name == "treatment_rep1") |>
    dplyr::pull(counts)
  expect_identical(test_counts_left, 9L)
})
##
### -----------------------------------------------------------------------###
