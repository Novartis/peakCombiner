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
reduced_colnames <- c(
  "chrom", "start", "end", "width", "strand", "name", "center", "score"
)
##
required_colnames <- c(
  "chrom", "start", "end", "name", "score", "strand",
  "center", "sample_name"
)
output_colnames <- c(
  "chrom", "start", "end", "width", "strand", "input_names"
)
##
data(syn_data_tibble)
test_data <- syn_data_tibble
input_colnames <- colnames(test_data)
##
test_data_prepared <- prepare_input_regions(
  data = test_data
)
test_data_center_expand <- center_expand_regions(
  data = test_data_prepared,
  center_by = "center_column",
  expand_by = NULL
)
test_data_filtered <- filter_regions(
  data = test_data_center_expand,
  exclude_by_blacklist = "hg38",
  include_by_chromosome_name = c("chr1", "chr10", "chr2", "chr42"),
  include_above_score_cutoff = NULL,
  include_top_n_scoring = NULL
)
test_data_disjoin_filter <- cr_disjoin_filter(
  data = test_data_filtered,
  found_in_samples = 2
)
test_data_reduce <- cr_reduce(
  data = test_data_disjoin_filter
)
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
  data <- test_data_reduce |>
    dplyr::mutate(chrom = as.character(chrom))
  ##
  expect_equal(length(colnames(data)), 8)
  expect_identical(names(data), reduced_colnames)
  expect_true(is.character(data$chrom))
  expect_true(is.numeric(data$start))
  expect_true(is.numeric(data$end))
  expect_true(is.character(data$name))
  expect_true(sum(stringr::str_detect(data$name, "|")) > 0)
  expect_true(sum(stringr::str_detect(data$name, "_")) > 0)
  ##
})
##
test_that("Input data frame has the expected structure", {
  ##
  data <- test_data_filtered |>
    dplyr::mutate(chrom = as.character(chrom))
  ##
  expect_equal(length(colnames(data)), 8)
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
### Test Output
### -----------------------------------------------------------------------###
##
test_that("Output data frame is correct", {
  ##
  data <- test_data_overlap |>
    dplyr::mutate(chrom = as.character(chrom))
  ##
  expect_setequal(colnames(data), reduced_colnames)
  expect_equal(ncol(data), 8)
  ##
  expect_identical(class(data)[2], "tbl")
  ##
  expect_true(is.character(data$chrom))
  expect_true(is.numeric(data$start))
  expect_true(is.numeric(data$end))
  expect_true(is.character(data$name))
  ##
  expect_identical(nrow(data), as.integer(41))
  expect_identical(data$start[1], 150L)
  expect_identical(sum(data$width), 31341L)
  ##
})
##
### -----------------------------------------------------------------------###
