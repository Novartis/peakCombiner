# Example
# test_that("multiplication works", {
#  expect_equal(2 * 2, 4)
# })
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
required_colnames <- c(
  "chrom", "start", "end", "name", "score", "strand", "center", "sample_name"
)
##
test_data <- peakCombiner::syn_data_tibble
input_colnames <- colnames(test_data)
##
test_data_prepared <- prepare_input_regions(
  data = test_data
)
test_data_center_expand <- center_expand_regions(
  data = test_data_prepared,
  center_by = "center_column",
  expand_by = NULL
)
test_data_filtered <- filter_regions(
  data = test_data_center_expand,
  exclude_by_blacklist = "hg38", # "hg38",
  include_by_chromosome_name = NULL,
  include_above_score_cutoff = NULL,
  include_top_n_scoring = NULL
)

##
test_data_disjoin_filter <- cr_disjoin_filter(
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
  expect_no_error(cr_disjoin_filter(
    data = test_data_filtered,
    found_in_samples = 3
  ))
  expect_error(cr_disjoin_filter(
    data = test_data_filtered,
    found_in_samples = 0
  ), "Arg")
  expect_error(cr_disjoin_filter(
    data = test_data_filtered,
    found_in_samples = NULL
  ), "Arg")
  expect_error(cr_disjoin_filter(
    data = test_data_filtered,
    found_in_samples = NA
  ), )
  expect_error(cr_disjoin_filter(
    data = test_data_filtered,
    found_in_samples = c(1, 2, 3)
  ), "'")
  expect_error(cr_disjoin_filter(
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
  expect_identical(nrow(data), as.integer(106))
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
