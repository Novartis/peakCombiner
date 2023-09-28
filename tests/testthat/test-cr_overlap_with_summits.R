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
input_colnames <- c("chr", "start", "end", "width", "strand", "input_names")
##
required_colnames <- c(
  "chr", "start", "end", "name", "score", "strand",
  "center", "sample_name"
)
output_colnames <- c("chr", "start", "end", "width", "strand", "input_names")
##
test_data <- readr::read_tsv("/da/ONC/BFx/research/muckema1/discovery_brd9/analysis/combpeaksr/lists/synthetic_genomic_regions.bed", show_col_types = FALSE)
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
test_data_filtered <- filter_regions(
  data = test_data_center_expand,
  filter_by_blacklist = "hg38", # "hg38",
  filter_by_chromosome_names = c("chr1", "chr10", "chr2", "chr42"),
  filter_by_significance = NULL,
  filter_by_top_enriched = NULL
) |> suppressWarnings()
test_data_disjoin_filter <- cr_disjoin_filter(data = test_data_filtered)
test_data_reduce <- cr_reduce(data = test_data_disjoin_filter)
##
test_data_overlap <- cr_overlap_with_summits(
  data = test_data_reduce,
  input = test_data_filtered
)
##
### -----------------------------------------------------------------------###
### Test input
### -----------------------------------------------------------------------###
##
test_that("Input data frame has the expected structure", {
  ##
  data <- test_data_reduce
  ##
  expect_equal(length(colnames(data)), 6)
  expect_identical(names(data), input_colnames)
  expect_true(is.character(data$chr))
  expect_true(is.numeric(data$start))
  expect_true(is.numeric(data$end))
  expect_true(is.character(data$input_names))
  expect_true(sum(str_detect(data$input_names, "|")) > 0)
  expect_true(sum(str_detect(data$input_names, ";")) > 0)
  ##
})
##
test_that("Input data frame has the expected structure", {
  ##
  data <- test_data_filtered
  ##
  expect_equal(length(colnames(data)), 8)
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
### Test Output
### -----------------------------------------------------------------------###
##
test_that("Output data frame is correct", {
  ##
  data <- test_data_overlap
  ##
  expect_setequal(colnames(data), output_colnames)
  expect_equal(ncol(data), 6)
  ##
  expect_identical(class(data)[2], "tbl")
  ##
  expect_true(is.character(data$chr))
  expect_true(is.numeric(data$start))
  expect_true(is.numeric(data$end))
  expect_true(is.character(data$input_names))
  ##
  expect_identical(nrow(data), as.integer(9))
  expect_identical(data$start[1], 150)
  expect_identical(sum(data$width), 6800L)
  ##
})
##
### -----------------------------------------------------------------------###
