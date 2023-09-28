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
test_data_prepared_filtered <- collapse_summits(data_prepared = test_data)
restult_colnames <- colnames(test_data_prepared_filtered)
##
### -----------------------------------------------------------------------###
### Test input
### -----------------------------------------------------------------------###
##
test_that("Test if function works with correct input", {
  expect_no_error(collapse_summits(data_prepared = test_data))
})
##
### -----------------------------------------------------------------------###
##
test_that("Input data has the right number of columns", {
  expect_equal(length(input_colnames), 8)
})
##
test_that("Column names of input data are identical with required once.", {
  expect_identical(names(test_data), required_colnames)
})
##
### -----------------------------------------------------------------------###
### Test output
### -----------------------------------------------------------------------###
##
test_that("Column names of output data are identical with required once.", {
  expect_setequal(colnames(test_data_prepared_filtered), required_colnames)
})
##
test_that("Output data has the right number of columns", {
  expect_equal(ncol(test_data_prepared_filtered), 8)
})
##
test_that("Output data has the right class.", {
  expect_identical(class(test_data_prepared_filtered)[2], "tbl")
})
##
test_that("Output data has the correct mean value for the column 'center'.", {
  expect_equal(mean(test_data_prepared_filtered$center), 1942.2535)
})
##
test_that("Output data has the correct number of rows.", {
  expect_identical(nrow(test_data_prepared_filtered), as.integer(71))
})
##
### -----------------------------------------------------------------------###