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
test_data_center_expand <- center_expand_regions(
  data = test_data_prepared,
  center_by = "summit",
  expand_by = NULL
)
##
input_colnames <- colnames(test_data_center_expand)
##
filter_by_significance <- 40
##
test_data_filtered <- filter_by_significance(
  data_filtered = test_data_center_expand,
  filter_by_significance = filter_by_significance
)
##
result_colnames <- colnames(test_data_filtered)
##
### -----------------------------------------------------------------------###
### Test input
### -----------------------------------------------------------------------###
##
test_that("Test if function works with correct input", {
  expect_no_error(filter_by_significance(
    data_filtered = test_data_center_expand,
    filter_by_significance = filter_by_significance
  ))
})
##
### -----------------------------------------------------------------------###
##
test_that("Input data frame has the expected structure", {
  data <- test_data_filtered
  ##
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
### -----------------------------------------------------------------------###
##
test_that("Required parameter 'filter_by_significance' has expected structure", {
  expect_no_error(filter_by_significance(
    data_filtered = test_data_filtered,
    filter_by_significance = NULL
  ))
  expect_no_error(filter_by_significance(
    data_filtered = test_data_filtered,
    filter_by_significance = 0
  ))
  ##
  expect_error(filter_by_significance(
    data_filtered = test_data_filtered,
    filter_by_significance = NA
  ))
  expect_error(filter_by_significance(
    data_filtered = test_data_filtered,
    filter_by_significance = "nonexisting"
  ))
  expect_error(filter_by_significance(
    data_filtered = test_data_filtered,
    filter_by_significance = c(1, 2, 3)
  ))
  ##
})
##
### -----------------------------------------------------------------------###
### Test Output
### -----------------------------------------------------------------------###
##
test_that("Output data frame is correct", {
  data <- test_data_filtered
  ##
  expect_setequal(colnames(data), required_colnames)
  expect_equal(ncol(data), 8)
  ##
  expect_identical(class(data)[2], "tbl")
  ##
  expect_true(is.character(data$chr))
  expect_true(is.numeric(data$start))
  expect_true(is.numeric(data$end))
  expect_true(is.character(data$name))
  expect_true(is.numeric(data$score))
  expect_true(is.character(data$strand))
  expect_true(is.numeric(data$center))
  expect_true(is.character(data$sample_name))
  ##
  expect_equal(mean(data$center), 2026.8657)
  expect_identical(nrow(data), as.integer(67))
  expect_identical(data$start[1], 250)
})
##
### -----------------------------------------------------------------------###
