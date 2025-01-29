##
### -----------------------------------------------------------------------###
### Prepare data for testing
### -----------------------------------------------------------------------###
##
set.seed(1234)
##
input_colnames_pre <- c(
  "chrom", "start", "end", "name", "score", "strand",
  "center", "sample_name"
)
input_colnames_post <- c(
  "chrom", "start", "end", "name", "score", "strand",
  "center", "sample_name", "input_names"
)
output_colnames_pre <- c(
  "chrom", "start", "end", "name", "score", "strand",
  "center", "sample_name"
)
output_colnames_post <- c(
  "chrom", "start", "end", "name", "score", "strand",
  "center", "sample_name", "input_names"
)
##
data(syn_data_bed, package = "peakCombiner")
test_data <- syn_data_bed
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
restult_colnames <- colnames(test_data_center_expand)
##
test_data_filtered <- peakCombiner::filter_regions(
  data = test_data_center_expand,
  exclude_by_blacklist = "hg38", # "hg38",
  include_by_chromosome_name = NULL,
  include_above_score_cutoff = NULL,
  include_top_n_scoring = NULL,
  show_messages = TRUE
)
##
test_data_combined <- peakCombiner::combine_regions(
  data = test_data_filtered,
  found_in_samples = 2,
  combined_center = "nearest",
  annotate_with_input_names = TRUE,
  combined_sample_name = "combined"
)
##
test_data_combined_ce <- peakCombiner::center_expand_regions(
  data = test_data_combined,
  center_by = "center_column",
  expand_by = NULL
)
### -----------------------------------------------------------------------###
### Test input
### -----------------------------------------------------------------------###

testthat::test_that("Test if function works with pre-combined input", {
  testthat::expect_no_error(peakCombiner::center_expand_regions(
    data = test_data_prepared,
    center_by = "center_column",
    expand_by = NULL
  ))
})

testthat::test_that("Test if function works with post-combined input", {
  testthat::expect_no_error(peakCombiner::center_expand_regions(
    data = test_data_combined,
    center_by = "center_column",
    expand_by = NULL
  ))
})

### -----------------------------------------------------------------------###

test_that("Required input data has the expected structure", {
  data <- test_data_prepared
  
  expect_equal(length(names(data)), 8)
  expect_identical(names(data), input_colnames_pre)
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

test_that("Required input data has the expected structure", {
  data <- test_data_combined
  
  expect_equal(length(names(data)), 9)
  expect_identical(names(data), input_colnames_post)
  expect_true(is.character(data$chrom))
  expect_true(is.numeric(data$start))
  expect_true(is.numeric(data$end))
  expect_true(is.character(data$name))
  expect_true(is.numeric(data$score))
  expect_true(is.character(data$strand))
  expect_true(is.numeric(data$center))
  expect_true(sum(stringr::str_detect(data$input_names, "|")) > 0)
})

### -----------------------------------------------------------------------###

test_that("Required paramter 'center_by' has the expected structure/value", {
  expect_no_error(peakCombiner::center_expand_regions(
    data = test_data_prepared,
    center_by = "center_Column",
    expand_by = NULL
  ))
  expect_error(peakCombiner::center_expand_regions(
    data = test_data_prepared,
    center_by = c("center_column", "calculated_value"),
    expand_by = NULL
  ), "center_by")
  expect_error(peakCombiner::center_expand_regions(
    data = test_data_prepared,
    center_by = "nonexisting",
    expand_by = NULL
  ), "center_by")
  expect_error(peakCombiner::center_expand_regions(
    data = test_data_prepared,
    center_by = NULL,
    expand_by = NULL
  ), "center_by")
  expect_error(peakCombiner::center_expand_regions(
    data = test_data_prepared,
    center_by = NA,
    expand_by = NULL
  ), "center_by")
})

### -----------------------------------------------------------------------###

testthat::test_that("Required paramter expand_by has the expected structure/value", {
  testthat::expect_no_error(peakCombiner::center_expand_regions(
    data = test_data_prepared,
    center_by = "center_column",
    expand_by = NULL
  ))
  testthat::expect_error(peakCombiner::center_expand_regions(
    data = test_data_prepared,
    center_by = "column_value",
    expand_by = NA
  ), )
  testthat::expect_error(peakCombiner::center_expand_regions(
    data = test_data_prepared,
    center_by = "column_value",
    expand_by = c(1, 2, 3)
  ), )
  testthat::expect_error(peakCombiner::center_expand_regions(
    data = test_data_prepared,
    center_by = "column_value",
    expand_by = "nonexisting"
  ), )
})

### -----------------------------------------------------------------------###
### Test Output
### -----------------------------------------------------------------------###

test_that("Output data frame is correct for pre-combined", {
  data <- test_data_center_expand
  
  expect_setequal(colnames(data), output_colnames_pre)
  expect_equal(ncol(data), 8)
  
  expect_identical(class(data)[2], "tbl")
  
  expect_true(is.character(data$chrom))
  expect_true(is.numeric(data$start))
  expect_true(is.numeric(data$end))
  expect_true(is.character(data$name))
  expect_true(is.numeric(data$score))
  expect_true(is.character(data$strand))
  expect_true(is.numeric(data$center))
  expect_true(is.character(data$sample_name))
  
  expect_equal(mean(data$center), 2495.6827)
  expect_identical(nrow(data), as.integer(52))
  expect_identical(data$start[1], 100.5)
})

test_that("Output data frame is correct for post-combined", {
  data <- test_data_combined_ce
  
  expect_setequal(colnames(data), output_colnames_post)
  expect_equal(ncol(data), 9)
  
  expect_identical(class(data)[2], "tbl")
  
  expect_true(is.character(data$chrom))
  expect_true(is.numeric(data$start))
  expect_true(is.numeric(data$end))
  expect_true(is.character(data$name))
  expect_true(is.numeric(data$score))
  expect_true(is.character(data$strand))
  expect_true(is.numeric(data$center))
  expect_equal(mean(data$center), 2770.45)
  expect_identical(nrow(data), as.integer(10))
  expect_identical(data$start[1], 200)
  expect_identical(data$end[1], 900)
  expect_identical(data$end[1], 900)
})

test_that("Output data frame is correct for data_prepared", {
  ##
  data <- test_data_prepared
  result <- peakCombiner::center_expand_regions(
    data = data,
    center_by = "center_column",
    expand_by = NULL
  )
  ##
  expect_no_error(peakCombiner::center_expand_regions(
    data = data,
    center_by = "center_column",
    expand_by = NULL
  ))
  ##
  expect_identical(nrow(result), 52L)
})
##
test_that("Output data frame is correct for data_center_expand", {
  ##
  data <- test_data_center_expand
  result <- peakCombiner::center_expand_regions(
    data = data,
    center_by = "center_column",
    expand_by = NULL
  )
  ##
  expect_no_error(peakCombiner::center_expand_regions(
    data = data,
    center_by = "center_column",
    expand_by = NULL
  ))
  ##
  expect_identical(nrow(result), 52L)
})
##
test_that("Output data frame is correct for data_filtered", {
  ##
  data <- test_data_filtered
  result <- peakCombiner::center_expand_regions(
    data = data,
    center_by = "center_column",
    expand_by = NULL
  )
  ##
  expect_no_error(peakCombiner::center_expand_regions(
    data = data,
    center_by = "center_column",
    expand_by = NULL
  ))
  ##
  expect_identical(nrow(result), 52L)
})
##
test_that("Output data frame is correct for data_combined", {
  ##
  data <- test_data_combined
  result <- peakCombiner::center_expand_regions(
    data = data,
    center_by = "center_column",
    expand_by = NULL
  )
  ##
  expect_no_error(peakCombiner::center_expand_regions(
    data = data,
    center_by = "center_column",
    expand_by = NULL
  ))
  ##
  expect_identical(nrow(result), 10L)
})
##
### -----------------------------------------------------------------------###
##
