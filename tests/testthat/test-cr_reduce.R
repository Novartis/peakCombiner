##
### -----------------------------------------------------------------------###
### Prepare data for testing
### -----------------------------------------------------------------------###
## tweak the prepare_input_regions() function and re-load it
#devtools::load_all()
##
### -----------------------------------------------------------------------###
### Prepare data for testing
### -----------------------------------------------------------------------###
##
set.seed(1234)
##
input_colnames <- c(
  "chrom", "start", "end", "width", "strand", "revmap", "sample_name",
  "ranking_comb_ref", "name", "center", "score", "rowname_disjoin"
)
##
output_colnames <- c("chr", "start", "end", "width", "strand", "input_names")
##
data(syn_data_tibble)
test_data <- syn_data_tibble
##
test_data_prepared <- prepare_input_regions(
  data = test_data
)
##
test_data_center_expand <- center_expand_regions(
  data = test_data_prepared,
  center_by = "center_column",
  expand_by = NULL
)
##
test_data_filtered <- filter_regions(
  data = test_data_center_expand,
  exclude_by_blacklist = "hg38",
  include_by_chromosome_name = c("chr1", "chr10", "chr2", "chr42"),
  include_above_score_cutoff = NULL,
  include_top_n_scoring = NULL
)
##
test_data_disjoin_filter <- cr_disjoin_filter(
  data = test_data_filtered,
  found_in_samples = 2
)
##
test_data_reduce <- cr_reduce(
  data = test_data_disjoin_filter
)
##
output_colnames <- colnames(
  test_data_reduce
)
##
### -----------------------------------------------------------------------###
### Test input
### -----------------------------------------------------------------------###
##
test_that("Input data frame has the expected structure", {
  data <- test_data_disjoin_filter |>
    dplyr::mutate(chrom = as.character(chrom))
  ##
  expect_equal(length(names(data)), 12)
  expect_identical(names(data), input_colnames)
  expect_true(is.character(data$chrom))
  expect_true(is.numeric(data$start))
  expect_true(is.numeric(data$end))
  expect_true(is.character(data$name))
  expect_true(sum(stringr::str_detect(data$name, "|")) > 0)
})
##
### -----------------------------------------------------------------------###
### Test Output
### -----------------------------------------------------------------------###
##
test_that("Output data frame is correct", {
  data <- test_data_reduce |>
    dplyr::mutate(chrom = as.character(chrom))
  ##
  expect_setequal(colnames(data), output_colnames)
  expect_equal(ncol(data), 8)
  ##
  expect_identical(class(data)[2], "tbl")
  ##
  expect_true(is.character(data$chrom))
  expect_true(is.numeric(data$start))
  expect_true(is.numeric(data$end))
  expect_true(is.character(data$name))
  ##
  expect_identical(nrow(data), 45L)
  expect_identical(data$start[1], 150L)
  expect_identical(round(sum(data$width),0), 31745)
  ##
})
##
### -----------------------------------------------------------------------###
