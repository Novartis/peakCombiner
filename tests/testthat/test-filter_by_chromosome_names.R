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
keep_chromosomes <- c("chr1", "chr10", "chr42")
##
test_data_filtered <- filter_by_chromosome_names(
  data_filtered = test_data_center_expand,
  filter_by_chromosome_names = keep_chromosomes
)
##
result_colnames <- colnames(test_data_filtered)
##
### -----------------------------------------------------------------------###
### Test input
### -----------------------------------------------------------------------###
##
test_that("Test if function works with correct input", {
  expect_no_error(filter_by_chromosome_names(
    data_filtered = test_data_center_expand,
    filter_by_chromosome_names = keep_chromosomes
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
test_that("Required parameter 'filter_by_chromosome_names' has expected structure", {
  expect_no_error(filter_by_chromosome_names(
    data_filtered = test_data_filtered,
    filter_by_chromosome_names = NULL
  )) |> suppressWarnings()
  expect_no_error(filter_by_chromosome_names(
    data_filtered = test_data_filtered,
    filter_by_chromosome_names = "chr1"
  )) |> suppressWarnings()
  expect_no_error(filter_by_chromosome_names(
    data_filtered = test_data_filtered,
    filter_by_chromosome_names = keep_chromosomes
  )) |> suppressWarnings()
  ##
  expect_error(filter_by_chromosome_names(
    data_filtered = test_data_filtered,
    filter_by_chromosome_names = NA
  ))
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
  expect_equal(mean(data$center), 2992.68293)
  expect_identical(nrow(data), as.integer(41))
  expect_identical(data$start[1], 250)
})
##
### -----------------------------------------------------------------------###
