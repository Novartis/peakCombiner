##
### -----------------------------------------------------------------------###
### Prepare data for testing
### -----------------------------------------------------------------------###
##
set.seed(1234)
##
reduced_colnames <- c(
  "chrom", "start", "end", "width", "strand", "name", "center", "score"
)
##
required_colnames <- c(
  "chrom", "start", "end", "name", "score", "strand",
  "center", "sample_name"
)
output_colnames <- c(
  "chrom", "start", "end", "width", "strand", "input_names"
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
  includeByChromosomeName = c("chr1", "chr10", "chr2", "chr42"),
  includeAboveScoreCutoff = NULL,
  includeTopNScoring = NULL,
  outputFormat = "tibble"
)
test_data_disjoin_filter <- peakCombiner:::crDisjoinFilter(
  data = test_data_filtered,
  foundInSamples = 2
)
test_data_reduce <- peakCombiner:::crReduce(
  data = test_data_disjoin_filter
)
##
test_data_overlap <- peakCombiner:::crOverlapWithSummits(
  data = test_data_reduce,
  input = test_data_filtered
)
##
### -----------------------------------------------------------------------###
### Test input
### -----------------------------------------------------------------------###
##
test_that("Input data frame has the expected structure", {
  ##
  data <- test_data_reduce |>
    dplyr::mutate(chrom = as.character(chrom))
  ##
  expect_equal(length(colnames(data)), 8)
  expect_identical(names(data), reduced_colnames)
  expect_true(is.character(data$chrom))
  expect_true(is.numeric(data$start))
  expect_true(is.numeric(data$end))
  expect_true(is.character(data$name))
  expect_true(sum(stringr::str_detect(data$name, "|")) > 0)
  expect_true(sum(stringr::str_detect(data$name, "_")) > 0)
  ##
})
##
test_that("Input data frame has the expected structure", {
  ##
  data <- test_data_filtered |>
    dplyr::mutate(chrom = as.character(chrom))
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
### -----------------------------------------------------------------------###
### Test Output
### -----------------------------------------------------------------------###
##
test_that("Output data frame is correct", {
  ##
  data <- test_data_overlap |>
    dplyr::mutate(chrom = as.character(chrom))
  ##
  expect_setequal(colnames(data), reduced_colnames)
  expect_equal(ncol(data), 8)
  ##
  expect_identical(class(data)[2], "tbl")
  ##
  expect_true(is.character(data$chrom))
  expect_true(is.numeric(data$start))
  expect_true(is.numeric(data$end))
  expect_true(is.character(data$name))
  ##
  expect_identical(nrow(data), as.integer(42))
  expect_identical(data$start[1], 1L)
  expect_identical(sum(data$width), 33300L)
  ##
})
##
### -----------------------------------------------------------------------###
