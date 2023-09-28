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
input_chr <- test_data_center_expand |>
  dplyr::pull(chr) |>
  unique()
##
test_keep_chromosomes <- extract_chromosome_names(
  data = test_data_center_expand,
  keep_chromosomes = "alphanumeric" # "alphanumeric"
)
##
### -----------------------------------------------------------------------###
### Test input
### -----------------------------------------------------------------------###
##
test_that("Test if function works with correct input", {
  expect_no_error(extract_chromosome_names(
    data = test_data_center_expand,
    keep_chromosomes = "alphanumeric"
  ))
})
##
### -----------------------------------------------------------------------###
##
test_that("Required colnumn names has the expected structure", {
  data <- test_data_center_expand
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
test_that("Required paramter 'data' has the expected structure/value", {
  expect_error(extract_chromosome_names(
    data = test_data_center_expand[2:4],
    keep_chromosomes = "all"
  ))
  expect_error(extract_chromosome_names(
    data = "nonexisting",
    keep_chromosomes = "all"
  ))
  expect_error(extract_chromosome_names(
    data = 1:10,
    keep_chromosomes = "all"
  ))
})
##
### -----------------------------------------------------------------------###
##
test_that("Required paramter 'center_by' has the expected structure/value", {
  expect_no_error(extract_chromosome_names(
    data = test_data_center_expand,
    keep_chromosomes = "aLphanumeric"
  ))
  expect_no_error(extract_chromosome_names(
    data = test_data_center_expand,
    keep_chromosomes = "ALL"
  ))
  ##
  ### -----------------------------------------------------------------------###
  ##
  expect_error(extract_chromosome_names(
    data = test_data_center_expand,
    keep_chromosomes = "chr1"
  ))
  expect_error(extract_chromosome_names(
    data = test_data_center_expand,
    keep_chromosomes = c("chr1", "chr10")
  ))
  expect_error(extract_chromosome_names(
    data = test_data_center_expand,
    keep_chromosomes = 1:10
  ))
  ##
})
##
### -----------------------------------------------------------------------###
### Test Output
### -----------------------------------------------------------------------###
##
test_that("Output data frame is correct", {
  expect_true(is.vector(test_keep_chromosomes))
  expect_true(all(test_keep_chromosomes %in% input_chr))
  ##
  expect_true(length(test_keep_chromosomes) == 3)
  expect_false("chr4 2" %in% test_keep_chromosomes)
})
##
### -----------------------------------------------------------------------###
