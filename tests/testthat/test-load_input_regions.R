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
input_data <- readr::read_tsv("/da/ONC/BFx/research/muckema1/discovery_brd9/analysis/opbaf-brd9-muckema1_rpackage_comb_peak/support/sample_sheet_test.tsv", show_col_types = FALSE)
# input_data <- readr::read_tsv("/da/ONC/BFx/research/muckema1/discovery_brd9/analysis/combpeaksr/lists/synthetic_genomic_regions.bed", show_col_types = FALSE)
input_colnames <- colnames(input_data)
##
all_colnames <- c(
  "chr", "start", "end", "name", "score", "strand", "center", "sample_name"
)
##
samplesheet_colnames <- c(
  "sample_name", "file_path", "file_format"
)
##
data_prepared <- load_input_regions(
  data = input_data,
  score_colname = "qValue",
  all_colnames = all_colnames
)
##
### -----------------------------------------------------------------------###
### Test input
### -----------------------------------------------------------------------###
##
test_that("Test if function works with correct input", {
  expect_no_error(load_input_regions(
    data = input_data,
    score_colname = NULL,
    all_colnames = all_colnames
  ))
})
##
test_that("Parameter `score_colname` has a length of 1.", {
  expect_equal(length("score_colname"), 1)
})
##
test_that("Parameter `score_colname` is a character vector.", {
  expect_equal(class("score_colname"), "character")
})
##
test_that("Parameter `score_colname` is not numeric.", {
  expect_false(is.numeric("score_colname"))
})
##
### -----------------------------------------------------------------------###
##
test_that("Input data has exact three columns.", {
  expect_equal(length(input_colnames), 3)
})
##
test_that("Input data colnames are the expected once.", {
  expect_identical(names(input_data), allowed_col_names)
})
##
test_that("Input column 'sample_name' is a class 'character'.", {
  expect_true(is.character(input_data$sample_name))
})
##
test_that("Input column 'file_path' is a class 'character'", {
  expect_true(is.character(input_data$file_path))
})
##
test_that("Input column 'file_format' is a class 'character'", {
  expect_true(is.character(input_data$file_format))
})
##
### -----------------------------------------------------------------------###
##
test_that("Error occurs when 'data' does not exist.", {
  expect_error(load_input_regions(
    data = "nonexisting",
    score_colname = "qValue",
    all_colnames = all_colnames
  ), "input_data")
})
##
test_that("Error occurs when 'data' has the wrong structure.", {
  expect_error(load_input_regions(
    data = tibble(1:10),
    score_colname = "qValue",
    all_colnames = all_colnames
  ), "input_data")
})
##
test_that("Error occurs when 'data' is a vector.", {
  expect_error(load_input_regions(
    data = as.vector(1:10),
    score_colname = "qValue",
    all_colnames = all_colnames
  ), "input_data")
})
##
test_that("Error occurs when 'data' is 'NULL'.", {
  expect_error(load_input_regions(
    data = NULL,
    score_colname = "qValue",
    all_colnames = all_colnames
  ), "input_data")
})
##
test_that("Error occurs when 'data' is 'NA'.", {
  expect_error(load_input_regions(
    data = NA,
    score_colname = "qValue",
    all_colnames = all_colnames
  ), "input_data")
})
##
### -----------------------------------------------------------------------###
##
test_that("Error occurs when 'score_colname' is not existing.", {
  expect_error(load_input_regions(
    data = input_data,
    score_colname = "nonexisting",
    all_colnames = all_colnames
  ), "score_colname")
})
##
test_that("Error occurs when 'score_colname' is other required colname.", {
  expect_error(load_input_regions(
    data = input_data,
    score_colname = "start",
    all_colnames = all_colnames
  ), "score_colname")
})
##
### -----------------------------------------------------------------------###
##
test_that("Error occurs when 'all_colnames' is not existing.", {
  expect_error(load_input_regions(
    input_data = input_data,
    score_colname = "qValue",
    all_colnames = "nonexisting"
  ), )
})
##
test_that("Error occurs when 'all_colnames' is vector of length 2.", {
  expect_error(load_input_regions(
    input_data = input_data,
    score_colname = "qValue",
    all_colnames = c("C1", "C2")
  ), )
})
##
test_that("Error occurs when 'all_colnames' is 'NULL'.", {
  expect_error(load_input_region(
    input_data = input_data,
    score_colname = "qValue",
    all_colnames = NULL
  ), )
})
##
### -----------------------------------------------------------------------###
### Test output
### -----------------------------------------------------------------------###
##
test_that("Column names of output data are identical with required once.", {
  expect_setequal(colnames(data_prepared), all_colnames)
})
##
test_that("Output data has the right number of columns", {
  expect_equal(ncol(data_prepared), 8)
  ##
})
##
test_that("Output data has the right class.", {
  expect_identical(class(data_prepared)[2], "tbl")
  ##
})
##
test_that("Output data has in column 'score', row 1 the correct value.", {
  expect_identical(data_prepared$score[1], 4701.96729)
})
##
test_that("Output data has the correct number of rows.", {
  expect_identical(nrow(data_prepared), 814153L)
})
