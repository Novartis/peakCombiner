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
  "chr", "start", "end", "name", "score", "strand",
  "center", "sample_name"
)
##
test_data <- readr::read_tsv("/da/ONC/BFx/research/muckema1/discovery_brd9/analysis/combpeaksr/lists/synthetic_genomic_regions.bed", show_col_types = FALSE)
input_colnames <- colnames(test_data)
##
test_data_prepared <- prepare_input_regions(
  input_data = test_data,
  score_colname = "qValue"
)
test_expansion_value <- define_expansion(
  data = test_data,
  expand_by = NULL
)
##
### -----------------------------------------------------------------------###
### Test input
### -----------------------------------------------------------------------###
##
test_that("Test if function works with correct input", {
  expect_no_error(define_expansion(
    data = test_data,
    expand_by = NULL
  ))
})
##
### -----------------------------------------------------------------------###
##
test_that("Required colnumn names has the expected structure", {
  data <- test_data
  expect_equal(length(input_colnames), 8)
  expect_identical(names(data), required_colnames)
  expect_true(is.character(data$chr))
  expect_true(is.numeric(data$start))
  expect_true(is.numeric(data$end))
  expect_true(is.character(data$name))
  expect_true(is.numeric(data$score))
  expect_true(is.character(data$strand))
  expect_true(is.numeric(data$center))
  expect_true(is.character(data$sample_name))
  expect_true(sum(str_detect(data$name, "|")) > 0)
})
##
test_that("Required colnumn names has the expected structure", {
  data <- test_data_prepared
  expect_equal(length(colnames(test_data_prepared)), 8)
  expect_identical(names(data), required_colnames)
  expect_true(is.character(data$chr))
  expect_true(is.numeric(data$start))
  expect_true(is.numeric(data$end))
  expect_true(is.character(data$name))
  expect_true(is.numeric(data$score))
  expect_true(is.character(data$strand))
  expect_true(is.numeric(data$center))
  expect_true(is.character(data$sample_name))
  expect_true(sum(str_detect(data$name, "|")) > 0)
})
##
### -----------------------------------------------------------------------###
##
test_that("Required paramter 'expand_by' has the expected structure", {
  expect_error(define_expansion(
    data = test_data,
    expand_by = NuLL
  ))
  expect_error(define_expansion(
    data = test_data,
    expand_by = "NULL"
  ))
  expect_error(define_expansion(
    data = test_data,
    expand_by = 0
  ))
  expect_error(define_expansion(
    data = test_data,
    expand_by = NA
  ))
  expect_error(define_expansion(
    data = test_data,
    expand_by = c(1, NA)
  ))
  expect_error(define_expansion(
    data = test_data,
    expand_by = c(1, "unexpected")
  ))
  expect_error(define_expansion(
    data = test_data,
    expand_by = c(1, 2, 3)
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
