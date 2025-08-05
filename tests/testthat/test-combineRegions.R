##
### -----------------------------------------------------------------------###
### Prepare data for testing
### -----------------------------------------------------------------------###
##
set.seed(1234)
##
input_colnames <- c(
  "chrom", "start", "end", "name", "score", "strand",
  "center", "sample_name"
)

output_colnames <- c(
  "chrom", "start", "end", "name", "score", "strand",
  "center", "sample_name", "center_origin", "input_names"
)

#' Prepare test data set
data(syn_data_tibble, package = "peakCombiner")
test_data <- syn_data_tibble
test_data

test_data_prepared <- peakCombiner:::prepareInputRegions(
  data = test_data,
  outputFormat = "tibble",
  showMessages = TRUE
)

test_data_center_expand <- peakCombiner:::centerExpandRegions(
  data = test_data_prepared,
  center_by = "center_column",
  expand_by = NULL,
  outputFormat = "tibble",
  showMessages = TRUE
)

test_data_filtered <- peakCombiner:::filterRegions(
  data = test_data_center_expand,
  include_by_chromosome_name = NULL,
  exclude_by_blacklist = NULL,
  include_above_score_cutoff = NULL,
  include_top_n_scoring = NULL,
  outputFormat = "tibble",
  showMessages = TRUE
)

test_data_combined <- peakCombiner:::combineRegions(
  data = test_data_filtered,
  combinedCenter = "nearest",
  annotateWithInputNames = FALSE,
  combinedSampleName = NULL,
  outputFormat = "tibble",
  showMessages = TRUE
)

### -----------------------------------------------------------------------###
### Test arguments
### -----------------------------------------------------------------------###

testthat::test_that("Input data frame has be data frame or tibble", {
  testthat::expect_error(peakCombiner:::combineRegions(
    data = c(1, 2, 3, 4, 5),
    outputFormat = "tibble",
    showMessages = FALSE
  ))
})

testthat::test_that("Input data frame has be data frame or tibble", {
  testthat::expect_error(peakCombiner:::combineRegions(
    data = NULL,
    outputFormat = "tibble",
    showMessages = FALSE
  ))
})

### -----------------------------------------------------------------------###
testthat::test_that("Argument 'combinedCenter' creates error if NULL", {
  testthat::expect_error(peakCombiner:::combineRegions(
    data = test_data_filtered,
    combinedCenter = NULL,
    outputFormat = "tibble",
    showMessages = FALSE
  ))
})

testthat::test_that("Argument 'combinedCenter' creates error if NA", {
  testthat::expect_error(peakCombiner:::combineRegions(
    data = test_data_filtered,
    combinedCenter = NA,
    outputFormat = "tibble",
    showMessages = FALSE
  ))
})

testthat::test_that("Argument 'combinedCenter' creates error if numeric
                    value", {
  testthat::expect_error(peakCombiner:::combineRegions(
    data = test_data_filtered,
    combinedCenter = 1,
    outputFormat = "tibble",
    showMessages = FALSE
  ))
})

testthat::test_that("Argument 'combinedCenter' tolerates capitilization", {
  testthat::expect_no_error(peakCombiner:::combineRegions(
    data = test_data_filtered,
    combinedCenter = "Nearest",
    outputFormat = "tibble",
    showMessages = FALSE
  ))
})

testthat::test_that("Argument 'combinedCenter' creates error if not allowes
                    value", {
  testthat::expect_error(peakCombiner:::combineRegions(
    data = test_data_filtered,
    combinedCenter = "Shortest",
    outputFormat = "tibble",
    showMessages = FALSE
  ))
})

### -----------------------------------------------------------------------###
testthat::test_that("Argument 'annotateWithInputNames' creates no error if
                    allowed value", {
  testthat::expect_no_error(peakCombiner:::combineRegions(
    data = test_data_filtered,
    annotateWithInputNames = TRUE,
    outputFormat = "tibble",
    showMessages = FALSE
  ))
  testthat::expect_no_error(peakCombiner:::combineRegions(
    data = test_data_filtered,
    annotateWithInputNames = FALSE,
    outputFormat = "tibble",
    showMessages = FALSE
  ))
})

testthat::test_that("Argument 'annotateWithInputNames' creates error if not
                    allowes value", {
  testthat::expect_error(peakCombiner:::combineRegions(
    data = test_data_filtered,
    annotateWithInputNames = FALSe,
    outputFormat = "tibble",
    showMessages = FALSE
  ))

  testthat::expect_error(peakCombiner:::combineRegions(
    data = test_data_filtered,
    annotateWithInputNames = 10,
    outputFormat = "tibble",
    showMessages = FALSE
  ))
})

testthat::test_that("Argument 'annotateWithInputNames' creates error if not
                    allowes value 'NA'", {
  testthat::expect_error(peakCombiner:::combineRegions(
    data = test_data_filtered,
    annotateWithInputNames = NA,
    outputFormat = "tibble",
    showMessages = FALSE
  ))
})

testthat::test_that("Argument 'annotateWithInputNames' creates error if not
                    allowes value 'NULL'", {
  testthat::expect_error(peakCombiner:::combineRegions(
    data = test_data_filtered,
    annotateWithInputNames = NULL,
    outputFormat = "tibble",
    showMessages = FALSE
  ))
})

testthat::test_that("Argument 'annotateWithInputNames' creates error if
                    length is greater then 1.", {
  testthat::expect_error(peakCombiner:::combineRegions(
    data = test_data_filtered,
    annotateWithInputNames = c(1, 2),
    outputFormat = "tibble",
    showMessages = FALSE
  ))
})

testthat::test_that("Argument 'annotateWithInputNames' creates error if not
                    allowed logical value with length 2 is provided.", {
  testthat::expect_error(peakCombiner:::combineRegions(
    data = test_data_filtered,
    annotateWithInputNames = c(NA, TRUE),
    outputFormat = "tibble",
    showMessages = FALSE
  ))
})
### -----------------------------------------------------------------------###

testthat::test_that("Argument 'combinedSampleName' creates no error if 'NULL'
          value is provided.", {
  testthat::expect_no_error(peakCombiner:::combineRegions(
    data = test_data_filtered,
    combinedSampleName = NULL,
    outputFormat = "tibble",
    showMessages = FALSE
  ))
})

testthat::test_that("Argument 'combinedSampleName' creates no error if single
                    character value is provided.", {
  testthat::expect_no_error(peakCombiner:::combineRegions(
    data = test_data_filtered,
    combinedSampleName = "Consensus",
    outputFormat = "tibble",
    showMessages = FALSE
  ))
})

testthat::test_that("Argument 'combinedSampleName' creates error if single
                    numeric value is provided.", {
  testthat::expect_error(peakCombiner:::combineRegions(
    data = test_data_filtered,
    combinedSampleName = 1,
    outputFormat = "tibble",
    showMessages = FALSE
  ))
})

testthat::test_that("Argument 'combinedSampleName' creates error if vector
                    with two entries is provided.", {
  testthat::expect_error(peakCombiner:::combineRegions(
    data = test_data_filtered,
    combinedSampleName = c("Consensus", "Two"),
    outputFormat = "tibble",
    showMessages = FALSE
  ))
})

testthat::test_that("Argument 'combinedSampleName' creates error if 'NA' is
          provided.", {
  testthat::expect_error(peakCombiner:::combineRegions(
    data = test_data_filtered,
    combinedSampleName = NA,
    outputFormat = "tibble",
    showMessages = FALSE
  ))
})

### -----------------------------------------------------------------------###

testthat::test_that("Argument 'showMessages' creates no error if TRUE or FALSE
          value is provided.", {
  testthat::expect_no_error(peakCombiner:::combineRegions(
    data = test_data_filtered,
    outputFormat = "tibble",
    showMessages = FALSE
  ))
  testthat::expect_no_error(peakCombiner:::combineRegions(
    data = test_data_filtered,
    outputFormat = "tibble",
    showMessages = TRUE
  ))
})

testthat::test_that("Argument 'showMessages' creates no error if non accepted
          value is provided.", {
  testthat::expect_error(peakCombiner:::combineRegions(
    data = test_data_filtered,
    outputFormat = "tibble",
    showMessages = FaLSE
  ))
})

testthat::test_that("Argument 'showMessages' creates no error if non accepted
          value 'NA' is provided.", {
  testthat::expect_error(peakCombiner:::combineRegions(
    data = test_data_filtered,
    outputFormat = "tibble",
    showMessages = NA
  ))
})


### -----------------------------------------------------------------------###
### Test input
### -----------------------------------------------------------------------###

testthat::test_that("Input data frame has the expected structure", {
  data <- test_data_filtered

  testthat::expect_equal(length(names(data)), 8)
  testthat::expect_identical(names(data), input_colnames)
  testthat::expect_true(is.character(data$chrom))
  testthat::expect_true(is.numeric(data$start))
  testthat::expect_true(is.numeric(data$end))
  testthat::expect_true(is.character(data$name))
  testthat::expect_true(is.numeric(data$score))
  testthat::expect_true(is.character(data$strand))
  testthat::expect_true(is.numeric(data$center))
  testthat::expect_true(is.character(data$sample_name))
  testthat::expect_true(sum(stringr::str_detect(data$name, "|")) > 0)
})

### -----------------------------------------------------------------------###
### Test Output
### -----------------------------------------------------------------------###

testthat::test_that("Output data has the correct classes and structure", {
  testthat::expect_no_error(peakCombiner:::check_data_structure(test_data_combined))
})

testthat::test_that("Output data frame has correct colnames", {
  testthat::expect_true(any(colnames(test_data_combined) %in% output_colnames))
})


### -----------------------------------------------------------------------###

testthat::test_that("Output data results has correct summit for 'nearest'
                    peak", {
  data <- peakCombiner:::combineRegions(
    data = test_data_filtered,
    foundInSamples = 2,
    combinedCenter = "nearest",
    annotateWithInputNames = FALSE,
    combinedSampleName = "consensus_peak",
    outputFormat = "tibble",
    showMessages = FALSE
  )

  testthat::expect_identical(round(data$center[7],0), 500)
  testthat::expect_identical(data$name[7], "consensus_peak|7")
})

test_that("Output data results has correct summit for 'strongst' peak", {
  data <- peakCombiner:::combineRegions(
    data = test_data_filtered,
    foundInSamples = 2,
    combinedCenter = "strongest",
    annotateWithInputNames = FALSE,
    combinedSampleName = "consensus_peak",
    outputFormat = "tibble",
    showMessages = FALSE
  )

  expect_identical(round(data$center[7],0), 600)
  expect_identical(data$name[7], "consensus_peak|7")
})

testthat::test_that("Output data results has correct summit for 'middle'
                    peak", {
  data <- peakCombiner:::combineRegions(
    data = test_data_filtered,
    foundInSamples = 2,
    combinedCenter = "middle",
    annotateWithInputNames = FALSE,
    combinedSampleName = "consensus_peak",
    outputFormat = "tibble",
    showMessages = FALSE
  )

  testthat::expect_identical(data$center[7], 550)
  testthat::expect_identical(data$name[7], "consensus_peak|7")
})
### -----------------------------------------------------------------------###
