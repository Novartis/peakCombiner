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
data(syn_data_tibble)
test_data <- syn_data_tibble
test_data

test_data_prepared <- prepare_input_regions(
  data = test_data,
  show_messages = TRUE
)

test_data_center_expand <- center_expand_regions(
  data = test_data_prepared,
  center_by = "center_column",
  expand_by = NULL,
  show_messages = TRUE
)

test_data_filtered <- filter_regions(
  data = test_data_center_expand,
  include_by_chromosome_name = NULL,
  exclude_by_blacklist = "hg38", # "hg38",
  include_above_score_cutoff = NULL,
  include_top_n_scoring = NULL,
  show_messages = TRUE
)

test_data_combined <- combine_regions(
  data = test_data_filtered,
  combined_center = "nearest",
  annotate_with_input_names = FALSE,
  combined_sample_name = NULL,
  show_messages = TRUE
)

### -----------------------------------------------------------------------###
### Test arguments
### -----------------------------------------------------------------------###

testthat::test_that("Input data frame has be data frame or tibble", {
  testthat::expect_error(combine_regions(
    data = c(1, 2, 3, 4, 5),
    show_messages = FALSE
  ))
})

testthat::test_that("Input data frame has be data frame or tibble", {
  testthat::expect_error(combine_regions(
    data = NULL,
    show_messages = FALSE
  ))
})

### -----------------------------------------------------------------------###
testthat::test_that("Argument 'combined_center' creates error if NULL", {
  testthat::expect_error(combine_regions(
    data = test_data_filtered,
    combined_center = NULL,
    show_messages = FALSE
  ))
})

testthat::test_that("Argument 'combined_center' creates error if NA", {
  testthat::expect_error(combine_regions(
    data = test_data_filtered,
    combined_center = NA,
    show_messages = FALSE
  ))
})

testthat::test_that("Argument 'combined_center' creates error if numeric
                    value", {
  testthat::expect_error(combine_regions(
    data = test_data_filtered,
    combined_center = 1,
    show_messages = FALSE
  ))
})

testthat::test_that("Argument 'combined_center' tolerates capitilization", {
  testthat::expect_no_error(combine_regions(
    data = test_data_filtered,
    combined_center = "Nearest",
    show_messages = FALSE
  ))
})

testthat::test_that("Argument 'combined_center' creates error if not allowes
                    value", {
  testthat::expect_error(combine_regions(
    data = test_data_filtered,
    combined_center = "Shortest",
    show_messages = FALSE
  ))
})

### -----------------------------------------------------------------------###
testthat::test_that("Argument 'annotate_with_input_names' creates no error if
                    allowed value", {
  testthat::expect_no_error(combine_regions(
    data = test_data_filtered,
    annotate_with_input_names = TRUE,
    show_messages = FALSE
  ))
  testthat::expect_no_error(combine_regions(
    data = test_data_filtered,
    annotate_with_input_names = FALSE,
    show_messages = FALSE
  ))
})

testthat::test_that("Argument 'annotate_with_input_names' creates error if not
                    allowes value", {
  testthat::expect_error(combine_regions(
    data = test_data_filtered,
    annotate_with_input_names = FALSe,
    show_messages = FALSE
  ))

  testthat::expect_error(combine_regions(
    data = test_data_filtered,
    annotate_with_input_names = 10,
    show_messages = FALSE
  ))
})

testthat::test_that("Argument 'annotate_with_input_names' creates error if not
                    allowes value 'NA'", {
  testthat::expect_error(combine_regions(
    data = test_data_filtered,
    annotate_with_input_names = NA,
    show_messages = FALSE
  ))
})

testthat::test_that("Argument 'annotate_with_input_names' creates error if not
                    allowes value 'NULL'", {
  testthat::expect_error(combine_regions(
    data = test_data_filtered,
    annotate_with_input_names = NULL,
    show_messages = FALSE
  ))
})

testthat::test_that("Argument 'annotate_with_input_names' creates error if
                    length is greater then 1.", {
  testthat::expect_error(combine_regions(
    data = test_data_filtered,
    annotate_with_input_names = c(1, 2),
    show_messages = FALSE
  ))
})

testthat::test_that("Argument 'annotate_with_input_names' creates error if not
                    allowed logical value with length 2 is provided.", {
  testthat::expect_error(combine_regions(
    data = test_data_filtered,
    annotate_with_input_names = c(NA, TRUE),
    show_messages = FALSE
  ))
})
### -----------------------------------------------------------------------###

testthat::test_that("Argument 'combined_sample_name' creates no error if 'NULL'
          value is provided.", {
  testthat::expect_no_error(combine_regions(
    data = test_data_filtered,
    combined_sample_name = NULL,
    show_messages = FALSE
  ))
})

testthat::test_that("Argument 'combined_sample_name' creates no error if single
                    character value is provided.", {
  testthat::expect_no_error(combine_regions(
    data = test_data_filtered,
    combined_sample_name = "Consensus",
    show_messages = FALSE
  ))
})

testthat::test_that("Argument 'combined_sample_name' creates error if single
                    numeric value is provided.", {
  testthat::expect_error(combine_regions(
    data = test_data_filtered,
    combined_sample_name = 1,
    show_messages = FALSE
  ))
})

testthat::test_that("Argument 'combined_sample_name' creates error if vector
                    with two entries is provided.", {
  testthat::expect_error(combine_regions(
    data = test_data_filtered,
    combined_sample_name = c("Consensus", "Two"),
    show_messages = FALSE
  ))
})

testthat::test_that("Argument 'combined_sample_name' creates error if 'NA' is
          provided.", {
  testthat::expect_error(combine_regions(
    data = test_data_filtered,
    combined_sample_name = NA,
    show_messages = FALSE
  ))
})

### -----------------------------------------------------------------------###

testthat::test_that("Argument 'show_messages' creates no error if TRUE or FALSE
          value is provided.", {
  testthat::expect_no_error(combine_regions(
    data = test_data_filtered,
    show_messages = FALSE
  ))
  testthat::expect_no_error(combine_regions(
    data = test_data_filtered,
    show_messages = TRUE
  ))
})

testthat::test_that("Argument 'show_messages' creates no error if non accepted
          value is provided.", {
  testthat::expect_error(combine_regions(
    data = test_data_filtered,
    show_messages = FaLSE
  ))
})

testthat::test_that("Argument 'show_messages' creates no error if non accepted
          value 'NA' is provided.", {
  testthat::expect_error(combine_regions(
    data = test_data_filtered,
    show_messages = NA
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
  testthat::expect_no_error(check_data_structure(test_data_combined))
})

testthat::test_that("Output data frame has correct colnames", {
  testthat::expect_true(any(colnames(test_data_combined) %in% output_colnames))
})


### -----------------------------------------------------------------------###

testthat::test_that("Output data results has correct summit for 'nearest'
                    peak", {
  data <- combine_regions(
    data = test_data_filtered,
    found_in_samples = 2,
    combined_center = "nearest",
    annotate_with_input_names = FALSE,
    combined_sample_name = "consensus_peak",
    show_messages = FALSE
  )

  testthat::expect_identical(round(data$center[7],0), 500)
  testthat::expect_identical(data$name[7], "consensus_peak|7")
})

test_that("Output data results has correct summit for 'strongst' peak", {
  data <- combine_regions(
    data = test_data_filtered,
    found_in_samples = 2,
    combined_center = "strongest",
    annotate_with_input_names = FALSE,
    combined_sample_name = "consensus_peak",
    show_messages = FALSE
  )

  expect_identical(round(data$center[7],0), 600)
  expect_identical(data$name[7], "consensus_peak|7")
})

testthat::test_that("Output data results has correct summit for 'middle'
                    peak", {
  data <- combine_regions(
    data = test_data_filtered,
    found_in_samples = 2,
    combined_center = "middle",
    annotate_with_input_names = FALSE,
    combined_sample_name = "consensus_peak",
    show_messages = FALSE
  )

  testthat::expect_identical(data$center[7], 550)
  testthat::expect_identical(data$name[7], "consensus_peak|7")
})
### -----------------------------------------------------------------------###
