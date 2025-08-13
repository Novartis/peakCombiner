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
test_data_prepared <- peakCombiner::prepareInputRegions(
  data = test_data,
  outputFormat = "tibble",
  showMessages = FALSE
)

test_data_center_expand <- peakCombiner::centerExpandRegions(
  data = test_data_prepared,
  centerBy = "center_column",
  outputFormat = "tibble",
  expandBy = NULL
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
  foundInSamples = 2)
test_data_reduce <- peakCombiner:::crReduce(data = test_data_disjoin_filter)
test_data_overlap <- peakCombiner:::crOverlapWithSummits(
  data = test_data_reduce,
  input = test_data_filtered
)
##
test_data_combined_with_summit <- peakCombiner:::crAddSummit(
  data = test_data_overlap,
  input = test_data_filtered,
  combinedCenter = "nearest",
  annotateWithInputNames = TRUE,
  combinedSampleName = "combined"
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
  expect_no_error(peakCombiner:::crAddSummit(
    data = test_data_overlap,
    input = test_data_filtered,
    combinedCenter = "STRONGEST"
  ))
  expect_no_error(peakCombiner:::crAddSummit(
    data = test_data_overlap,
    input = test_data_filtered,
    combinedCenter = "middle"
  ))
  ##
  expect_error(peakCombiner:::crAddSummit(
    data = test_data_overlap,
    input = test_data_filtered,
    combinedCenter = mean
  ))
  expect_error(peakCombiner:::crAddSummit(
    data = test_data_overlap,
    input = test_data_filtered,
    combinedCenter = 2
  ))
  expect_error(peakCombiner:::crAddSummit(
    data = test_data_overlap,
    input = test_data_filtered,
    combinedCenter = c(1, 2, 3)
  ), "`")
  expect_error(peakCombiner:::crAddSummit(
    data = test_data_overlap,
    input = test_data_filtered,
    combinedCenter = NULL
  ), "`")
  expect_error(peakCombiner:::crAddSummit(
    data = test_data_overlap,
    input = test_data_filtered,
    combinedCenter = NA
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
  expect_identical(nrow(data), as.integer(3))
  expect_identical(data$center[1], 301)
  expect_identical(round(sum(data$score), 0), 245)
  ##
})
##
test_that("Output data results with different summits", {
  data <- peakCombiner:::crAddSummit(
    data = test_data_overlap,
    input = test_data_filtered,
    combinedCenter = "nearest"
  )
  expect_identical(data$center[2], 301)
  ##
  data <- peakCombiner:::crAddSummit(
    data = test_data_overlap,
    input = test_data_filtered,
    combinedCenter = "strongest"
  )
  expect_identical(data$center[2], 301)
  ##
  data <- peakCombiner:::crAddSummit(
    data = test_data_overlap,
    input = test_data_filtered,
    combinedCenter = "middle"
  )
  expect_identical(data$center[2], 325.5)
  ##
})
##
