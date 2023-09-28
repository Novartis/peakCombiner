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
  "chr", "start", "end", "width", "strand", "revmap",
  "ranking_comb_ref", "name", "rowname_disjoin"
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
test_data_filtered <- filter_regions(
  data = test_data_center_expand,
  filter_by_blacklist = "hg38", # "hg38",
  filter_by_chromosome_names = NULL,
  filter_by_significance = NULL,
  filter_by_top_enriched = NULL
) |> suppressWarnings()
##
input_colnames <- colnames(test_data_filtered)
##
test_data_disjoin_filter <- cr_disjoin_filter(
  data = test_data_filtered,
  found_in_samples = 2
)
##
result_colnames <- colnames(test_data_disjoin_filter)
##
### -----------------------------------------------------------------------###
### Test input
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
test_that("Parameter 'found_in_samples' has the correct structure", {
  expect_no_error(cr_disjoin_filter(
    data = test_data_filtered,
    found_in_samples = 3
  ))
  expect_error(cr_disjoin_filter(
    data = test_data_filtered,
    found_in_samples = 0
  ), "'")
  expect_error(cr_disjoin_filter(
    data = test_data_filtered,
    found_in_samples = NULL
  ), "'")
  expect_error(cr_disjoin_filter(
    data = test_data_filtered,
    found_in_samples = NA
  ), "'")
  expect_error(cr_disjoin_filter(
    data = test_data_filtered,
    found_in_samples = c(1, 2, 3)
  ), "'")
  expect_error(cr_disjoin_filter(
    data = test_data_filtered,
    found_in_samples = test_data_filtered
  ), "'")
})
### -----------------------------------------------------------------------###
### Test Output
### -----------------------------------------------------------------------###
##
test_that("Output data frame is correct", {
  data <- test_data_disjoin_filter
  ##
  expect_setequal(colnames(data), required_colnames)
  expect_equal(ncol(data), 9)
  ##
  expect_identical(class(data)[2], "tbl")
  ##
  expect_true(is.character(data$chr))
  expect_true(is.numeric(data$start))
  expect_true(is.numeric(data$end))
  expect_true(is.character(data$sample_name))
  ##
  expect_identical(nrow(data), as.integer(186))
  expect_identical(data$start[1], 150)
  ##
  test_counts_left <- test_data_filtered |>
    dplyr::group_by(sample_name) |>
    dplyr::summarise(counts = n()) |>
    dplyr::filter(sample_name == "treatment_rep1") |>
    dplyr::pull(counts)
  expect_identical(test_counts_left, as.integer(12))
})
##
### -----------------------------------------------------------------------###
