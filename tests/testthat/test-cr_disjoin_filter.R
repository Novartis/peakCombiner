##
### -----------------------------------------------------------------------###
### Prepare data for testing
### -----------------------------------------------------------------------###
##
set.seed(1234)
##
required_colnames <- c(
  "chrom", "start", "end", "name", "score", "strand", "center", "sample_name"
)
##
data(syn_data_tibble, package = "peakCombiner")
test_data <- syn_data_tibble
input_colnames <- colnames(test_data)
##
test_data_prepared <- peakCombiner::prepareInputRegions(
  data = test_data,
  outputFormat = "tibble",
  showMessages = FALSE
)
test_data_center_expand <- peakCombiner::centerExpandRegions(
  data = test_data_prepared,
  centerBy = "center_column",
  expandBy = NULL,
  outputFormat = "tibble"
)
test_data_filtered <- peakCombiner::filterRegions(
  data = test_data_center_expand,
  excludeByBlacklist = NULL,
  includeByChromosomeName = NULL,
  includeAboveScoreCutoff = NULL,
  includeTopNScoring = NULL,
  outputFormat = "tibble"
)

##
test_data_disjoin_filter <- peakCombiner:::crDisjoinFilter(
  data = test_data_filtered,
  foundInSamples = 2
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
test_that("Parameter 'foundInSamples' has the correct structure", {
  expect_no_error(peakCombiner:::crDisjoinFilter(
    data = test_data_filtered,
    foundInSamples = 3
  ))
  expect_error(peakCombiner:::crDisjoinFilter(
    data = test_data_filtered,
    foundInSamples = 0
  ), "Arg")
  expect_error(peakCombiner:::crDisjoinFilter(
    data = test_data_filtered,
    foundInSamples = NULL
  ), "Arg")
  expect_error(peakCombiner:::crDisjoinFilter(
    data = test_data_filtered,
    foundInSamples = NA
  ), )
  expect_error(peakCombiner:::crDisjoinFilter(
    data = test_data_filtered,
    foundInSamples = c(1, 2, 3)
  ), "'")
  expect_error(peakCombiner:::crDisjoinFilter(
    data = test_data_filtered,
    foundInSamples = test_data_filtered
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
  expect_identical(nrow(data), as.integer(162))
  expect_identical(data$start[1], 1)
  ##
  test_counts_left <- test_data_filtered |>
    dplyr::group_by(sample_name) |>
    dplyr::summarise(counts = dplyr::n()) |>
    dplyr::filter(sample_name == "treatment_rep1") |>
    dplyr::pull(counts)
  expect_identical(test_counts_left, as.integer(8))
})
##
### -----------------------------------------------------------------------###
