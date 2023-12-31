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
input_colnames_pre <- c(
  "chrom", "start", "end", "name", "score", "strand",
  "center", "sample_name"
)
input_colnames_post <- c(
  "chrom", "start", "end", "name", "score", "strand",
  "center", "sample_name", "center_origin", "input_names"
)
output_colnames_pre <- c(
  "chrom", "start", "end", "name", "score", "strand",
  "center", "sample_name"
)
output_colnames_post <- c(
  "chrom", "start", "end", "name", "score", "strand",
  "center", "sample_name", "center_origin", "input_names"
)
##
test_data <- readr::read_tsv("lists/synthetic_genomic_regions.bed", show_col_types = FALSE)
##
test_data_prepared <- prepare_input_regions(
  data = test_data
  )
##
test_data_center_expand <- center_expand_regions(
  data = test_data_prepared,
  center_by = "column_value",
  expand_by = NULL
)
restult_colnames <- colnames(test_data_center_expand)
##
test_data_filtered <- filter_regions(
  data = test_data_center_expand,
  filter_by_blacklist = "hg38", # "hg38",
  filter_by_chromosome_names = NULL,
  filter_by_significance = NULL,
  filter_by_top_enriched = NULL
) |> suppressWarnings()
##
test_data_combined <- combine_regions(
  data = test_data_filtered,
  found_in_samples = 2,
  center = "nearest"
)
##
test_data_combined_ce <- center_expand_regions(
  data = test_data_combined,
  center_by = "column_value",
  expand_by = NULL
)
### -----------------------------------------------------------------------###
### Test input
### -----------------------------------------------------------------------###

test_that("Test if function works with pre-combined input", {
  expect_no_error(center_expand_regions(
    data = test_data_prepared,
    center_by = "column_value",
    expand_by = NULL
  ))
})

test_that("Test if function works with post-combined input", {
  expect_no_error(center_expand_regions(
    data = test_data_combined,
    center_by = "column_value",
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
  expect_true(sum(str_detect(data$name, "|")) > 0)
})

test_that("Required input data has the expected structure", {
  
  data <- test_data_combined
  
  expect_equal(length(names(data)), 10)
  expect_identical(names(data), input_colnames_post)
  expect_true(is.character(data$chrom))
  expect_true(is.numeric(data$start))
  expect_true(is.numeric(data$end))
  expect_true(is.character(data$name))
  expect_true(is.numeric(data$score))
  expect_true(is.character(data$strand))
  expect_true(is.numeric(data$center))
  expect_true(is.character(data$center_origin))
  expect_true(sum(str_detect(data$input_names, "|")) > 0)
})

### -----------------------------------------------------------------------###

test_that("Required paramter 'center_by' has the expected structure/value", {
  expect_no_error(center_expand_regions(
    data = test_data_prepared,
    center_by = "coluMn_value",
    expand_by = NULL
  ))
  expect_error(center_expand_regions(
    data = test_data_prepared,
    center_by = c("column_value", "calculated_value"),
    expand_by = NULL
  ), "center_by")
  expect_error(center_expand_regions(
    data = test_data_prepared,
    center_by = "nonexisting",
    expand_by = NULL
  ), "center_by")
  expect_error(center_expand_regions(
    data = test_data_prepared,
    center_by = NULL,
    expand_by = NULL
  ), "center_by")
  expect_error(center_expand_regions(
    data = test_data_prepared,
    center_by = NA,
    expand_by = NULL
  ), "center_by")
})

### -----------------------------------------------------------------------###

test_that("Required paramter expand_by has the expected structure/value", {
  expect_no_error(center_expand_regions(
    data = test_data_prepared,
    center_by = "column_value",
    expand_by = NULL
  ))
  expect_error(center_expand_regions(
    data = test_data_prepared,
    center_by = "column_value",
    expand_by = NA
  ), )
  expect_error(center_expand_regions(
    data = test_data_prepared,
    center_by = "column_value",
    expand_by = c(1, 2, 3)
  ), )
  expect_error(center_expand_regions(
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

  expect_equal(mean(data$center), 1942.2535)
  expect_identical(nrow(data), as.integer(71))
  expect_identical(data$start[1], 250)
})

test_that("Output data frame is correct for post-combined", {

  data <- test_data_combined_ce

  expect_setequal(colnames(data), output_colnames_post)
  expect_equal(ncol(data), 10)

  expect_identical(class(data)[2], "tbl")

  expect_true(is.character(data$chrom))
  expect_true(is.numeric(data$start))
  expect_true(is.numeric(data$end))
  expect_true(is.character(data$name))
  expect_true(is.numeric(data$score))
  expect_true(is.character(data$strand))
  expect_true(is.numeric(data$center))
  expect_true(is.character(data$center_origin))

  expect_equal(mean(data$center), 2192.3077)
  expect_identical(nrow(data), as.integer(13))
  expect_identical(data$start[1], 100)
  expect_identical(data$end[1], 900)
  expect_identical(data$end[1], 900)
})

test_that("Output data frame is correct for data_prepared", {
  ##
  data <- test_data_prepared
  result <- center_expand_regions(
    data = data,
    center_by = "column_value",
    expand_by = NULL
  )
  ##
  expect_no_error(center_expand_regions(
    data = data,
    center_by = "column_value",
    expand_by = NULL
  ))
  ##
  expect_identical(nrow(result), 71L)
  expect_identical(result$start[9], as.numeric(250))
})
##
test_that("Output data frame is correct for data_center_expand", {
  ##
  data <- test_data_center_expand
  result <- center_expand_regions(
    data = data,
    center_by = "column_value",
    expand_by = NULL
  )
  ##
  expect_no_error(center_expand_regions(
    data = data,
    center_by = "column_value",
    expand_by = NULL
  ))
  ##
  expect_identical(nrow(result), 71L)
  expect_identical(result$start[9], as.numeric(250))
})
##
test_that("Output data frame is correct for data_filtered", {
  ##
  data <- test_data_filtered
  result <- center_expand_regions(
    data = data,
    center_by = "column_value",
    expand_by = NULL
  )
  ##
  expect_no_error(center_expand_regions(
    data = data,
    center_by = "column_value",
    expand_by = NULL
  ))
  ##
  expect_identical(nrow(result), 71L)
  expect_identical(result$start[9], as.numeric(250))
})
##
test_that("Output data frame is correct for data_combined", {
  ##
  data <- test_data_combined
  result <- center_expand_regions(
    data = data,
    center_by = "column_value",
    expand_by = NULL
  )
  ##
  expect_no_error(center_expand_regions(
    data = data,
    center_by = "column_value",
    expand_by = NULL
  ))
  ##
  expect_identical(nrow(result), 13L)
  expect_identical(result$start[9], as.numeric(100))
})
##
### -----------------------------------------------------------------------###
##
test_that("All newly center and expand regions do have only positive values", {
  ##
  data_input <- readr::read_tsv(paste0("lists/input_data-bed3.bed"))
  prep_data <- prepare_input_regions(data = data_input,score_colname = NULL)
  ##
  expect_no_error(center_expand_regions(
    data = prep_data,
    center_by = "column_value",
    expand_by = NULL
  ))
  ##
  expect_message(center_expand_regions(
    data = prep_data,
    center_by = "column_value",
    expand_by = 500
  ))
})
