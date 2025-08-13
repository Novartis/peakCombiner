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
test_expansion_value <- 350
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
##
### -----------------------------------------------------------------------###
### Test input
### -----------------------------------------------------------------------###
##
test_that("Test if function works with correct input", {
  expect_no_error(peakCombiner:::defineExpansion(
    data = test_data,
    expandBy = NULL
  ))
})
##
### -----------------------------------------------------------------------###
##
test_that("Required colnumn names has the expected structure", {
  data <- test_data

  expect_equal(length(input_colnames), 7)
  expect_identical(names(test_data_prepared), required_colnames)
  expect_true(is.character(data$chrom))
  expect_true(is.numeric(data$start))
  expect_true(is.numeric(data$end))
  expect_true(is.numeric(data$score))
  expect_true(is.character(data$strand))
  #expect_true(is.numeric(data$center))
  # expect_true(sum(stringr::str_detect(data$name, "|")) > 0)
})
##
test_that("Required colnumn names has the expected structure", {
  data <- test_data_prepared
  expect_equal(length(colnames(test_data_prepared)), 8)
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
test_that("Required paramter 'expand_by' has the expected structure", {
  expect_error(define_expansion(
    data = test_data,
    expandBy = NuLL
  ))
  expect_error(define_expansion(
    data = test_data,
    expandBy = "NULL"
  ))
  expect_error(define_expansion(
    data = test_data,
    expandBy = 0
  ))
  expect_error(define_expansion(
    data = test_data,
    expandBy = NA
  ))
  expect_error(define_expansion(
    data = test_data,
    expandBy = c(1, NA)
  ))
  expect_error(define_expansion(
    data = test_data,
    expandBy = c(1, "unexpected")
  ))
  expect_error(define_expansion(
    data = test_data,
    expandBy = c(1, 2, 3)
  ))
})
##
### -----------------------------------------------------------------------###
### Test Output
### -----------------------------------------------------------------------###
##
test_that("Output has the expected structure", {
  expect_true(is.numeric(test_expansion_value))
  expect_true(test_expansion_value > 0)
  expect_false(length(test_expansion_value) > 2)
  expect_identical(test_expansion_value, 350)
})
##
