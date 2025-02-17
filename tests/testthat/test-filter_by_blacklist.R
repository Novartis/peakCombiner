##
### -----------------------------------------------------------------------###
### Prepare data for testing
### -----------------------------------------------------------------------###
##
set.seed(1234)
##
required_colnames <- c(
  "chrom", "start", "end", "name", "score", "strand",
  "center", "sample_name"
)
##
data(syn_data_tibble, package = "peakCombiner")
test_data <- syn_data_tibble
input_colnames <- colnames(test_data)
##
test_data_prepared <- peakCombiner::prepare_input_regions(
  data = test_data
)
##
test_data_center_expand <- peakCombiner::center_expand_regions(
  data = test_data_prepared,
  center_by = "center_column",
  expand_by = NULL
)
##
test_data_filtered <- peakCombiner:::filter_by_chromosome_names(
  data = test_data_center_expand,
  include_by_chromosome_name = c("chr1", "chr10", "chr42")
)
##
input_colnames <- colnames(test_data_filtered)
##
data(blacklist_hg38, package = "peakCombiner")
blacklist <- blacklist_hg38
##
test_data_filtered_bl <- peakCombiner:::filter_by_blacklist(
  data = test_data_filtered,
  exclude_by_blacklist = blacklist
)
##
result_colnames <- colnames(test_data_filtered)
##
### -----------------------------------------------------------------------###
### Test input
### -----------------------------------------------------------------------###
##
##
test_that("Test if function works with correct input", {
  expect_no_error(peakCombiner:::filter_by_blacklist(
    data = test_data_filtered,
    exclude_by_blacklist = blacklist
  ))
  expect_no_error(peakCombiner:::filter_by_blacklist(
    data = test_data_filtered,
    exclude_by_blacklist = NULL
  ))
  expect_no_error(peakCombiner:::filter_by_blacklist(
    data = test_data_filtered,
    exclude_by_blacklist = "hg38"
  ))
})
##
### -----------------------------------------------------------------------###
##
test_that("Input data frame has the expected structure", {
  data <- test_data_filtered

  expect_equal(length(input_colnames), 8)
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
##
test_that("Required parameter 'filter_by_blacklist' has expected structure", {
  expect_no_error(peakCombiner:::filter_by_blacklist(
    data = test_data_filtered,
    exclude_by_blacklist = NULL
  ))
  expect_no_error(peakCombiner:::filter_by_blacklist(
    data = test_data_filtered,
    exclude_by_blacklist = "HG38"
  ))
  expect_no_error(peakCombiner:::filter_by_blacklist(
    data = test_data_filtered,
    exclude_by_blacklist = "mm10"
  ))
  ##
  expect_error(peakCombiner:::filter_by_blacklist(
    data = test_data_filtered,
    filter_by_blacklist = blacklist[1:2]
  ))
})
##
### -----------------------------------------------------------------------###
##
test_that("For 'filter_by_blacklist' providing blacklist with different
          names", {
  blacklist2 <- blacklist
  colnames(blacklist2) <- c("CHROM", "start", "end")
  ##
  expect_no_error(peakCombiner:::filter_by_blacklist(
    data = test_data_filtered,
    exclude_by_blacklist = blacklist2
  ))
  ##
  colnames(blacklist2) <- c("seqnames", "start", "end")
  ##
  expect_error(peakCombiner:::filter_by_blacklist(
    data = test_data_filtered,
    exclude_by_blacklist = blacklist2
  ))
})
##
### -----------------------------------------------------------------------###
##
test_that("Wrong input for exclude_by_blacklist for 'filter_by_blacklist'", {
  expect_error(peakCombiner:::filter_by_blacklist(
    data = test_data_filtered,
    exclude_by_blacklist = "mm38"
  ))
  expect_error(peakCombiner:::filter_by_blacklist(
    data = test_data_filtered,
    exclude_by_blacklist = hg38
  ))
  expect_error(peakCombiner:::filter_by_blacklist(
    data = test_data_filtered,
    exclude_by_blacklist = 1
  ))
  expect_error(peakCombiner:::filter_by_blacklist(
    data = test_data_filtered,
    exclude_by_blacklist = c(1, 2)
  ))
})
##
### -----------------------------------------------------------------------###
### Test Output
### -----------------------------------------------------------------------###
##
test_that("Output data frame is correct", {
  data <- test_data_filtered_bl
  ##
  expect_setequal(colnames(data), required_colnames)
  expect_equal(ncol(data), 8)
  ##
  expect_identical(class(data)[2], "tbl")
  ##
  expect_true(is.character(data$chrom))
  expect_true(is.numeric(data$start))
  expect_true(is.numeric(data$end))
  expect_true(is.character(data$name))
  expect_true(is.numeric(data$score))
  expect_true(is.character(data$strand))
  expect_true(is.numeric(data$center))
  expect_true(is.character(data$sample_name))
  ##
  expect_equal(round(mean(data$center), 0), 3168)
  expect_identical(nrow(data), 38L)
  expect_identical(data$start[1], 250L)
})
##
### -----------------------------------------------------------------------###
