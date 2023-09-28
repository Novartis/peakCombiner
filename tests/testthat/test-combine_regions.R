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

### -----------------------------------------------------------------------###
### Prepare data for testing
### -----------------------------------------------------------------------###

input_colnames <- c(
  "chrom", "start", "end", "name", "score", "strand",
  "center", "sample_name"
)

output_colnames <- c(
  "chrom", "start", "end", "name", "score", "strand",
  "center", "sample_name", "center_origin", "input_names"
)

test_data  <- readr::read_tsv(paste0("lists/synthetic_genomic_regions.bed"), show_col_types = FALSE)

test_data_prepared <- prepare_input_regions(
  data = test_data,
  show_messages = TRUE
)

test_data_center_expand <- center_expand_regions(
  data = test_data_prepared,
  center_by = "column_value",
  expand_by = NULL,
  show_messages = TRUE
)

test_data_filtered <- filter_regions(
  data = test_data_center_expand,
  include_by_chromosome_name = NULL,
  exclude_by_blacklist = "hg38", # "hg38",
  include_above_score_cutoff = NULL,
  include_top_n_scoring = NULL,
  show_messages = TRUE
)

test_data_combined <- combine_regions(
  data = test_data_filtered,
  combined_center = "nearest",
  annotate_with_input_names = FALSE,
  combined_sample_name = NULL,
  show_messages = TRUE
)

### -----------------------------------------------------------------------###
### Test arguments
### -----------------------------------------------------------------------###

test_that("Input data frame has be data frame or tibble", {
  expect_error(combine_regions(data = c(1,2,3,4,5),
                               show_messages = FALSE))
})

test_that("Input data frame has be data frame or tibble", {
  expect_error(combine_regions(data = NULL,
                               show_messages = FALSE))
})

### -----------------------------------------------------------------------###
test_that("Argument 'combined_center' creates error if NULL", {
  expect_error(combine_regions(data = test_data_filtered,
                               combined_center = NULL,
                               show_messages = FALSE))
})

test_that("Argument 'combined_center' creates error if NA", {
  expect_error(combine_regions(data = test_data_filtered,
                               combined_center = NA,
                               show_messages = FALSE))
})

test_that("Argument 'combined_center' creates error if numeric value", {
  expect_error(combine_regions(data = test_data_filtered,
                               combined_center = 1,
                               show_messages = FALSE))
})

test_that("Argument 'combined_center' tolerates capitilization", {
  expect_no_error(combine_regions(data = test_data_filtered,
                               combined_center = "Nearest",
                               show_messages = FALSE))
})

test_that("Argument 'combined_center' creates error if not allowes value", {
  expect_error(combine_regions(data = test_data_filtered,
                                  combined_center = "Shortest",
                                  show_messages = FALSE))
})

### -----------------------------------------------------------------------###
test_that("Argument 'annotate_with_input_names' creates no error if allowed 
          value", {
            expect_error(combine_regions(data = test_data_filtered,
                                         annotate_with_input_names = TRUE,
                                         show_messages = FALSE))
            expect_error(combine_regions(data = test_data_filtered,
                                         annotate_with_input_names = FALSE,
                                         show_messages = FALSE))
          })

test_that("Argument 'annotate_with_input_names' creates error if not allowes 
          value", {
  expect_error(combine_regions(data = test_data_filtered,
                               annotate_with_input_names = FALSe,
                               show_messages = FALSE))
})

test_that("Argument 'annotate_with_input_names' creates error if not allowes 
          value 'NA'", {
            expect_error(combine_regions(data = test_data_filtered,
                                         annotate_with_input_names = NA,
                                         show_messages = FALSE))
          })

test_that("Argument 'annotate_with_input_names' creates error if not allowes 
          value 'NULL'", {
            expect_error(combine_regions(data = test_data_filtered,
                                         annotate_with_input_names = NULL,
                                         show_messages = FALSE))
          })

test_that("Argument 'annotate_with_input_names' creates error if length is 
          greater then 1.", {
            expect_error(combine_regions(data = test_data_filtered,
                                         annotate_with_input_names = c(1,2),
                                         show_messages = FALSE))
          })

test_that("Argument 'annotate_with_input_names' creates error if not allowed 
logical value with length 2 is provided.", {
            expect_error(combine_regions(data = test_data_filtered,
                                         annotate_with_input_names = c(NA,TRUE),
                                         show_messages = FALSE))
          })
### -----------------------------------------------------------------------###

test_that("Argument 'combined_sample_name' creates no error if 'NULL' 
          value is provided.", {
            expect_no_error(combine_regions(data = test_data_filtered,
                                            combined_sample_name = NULL,
                                            show_messages = FALSE))
          })

test_that("Argument 'combined_sample_name' creates no error if single character 
          value is provided.", {
  expect_no_error(combine_regions(data = test_data_filtered,
                               combined_sample_name = "Consensus",
                               show_messages = FALSE))
})

test_that("Argument 'combined_sample_name' creates error if single numeric 
          value is provided.", {
            expect_error(combine_regions(data = test_data_filtered,
                                            combined_sample_name = 1,
                                            show_messages = FALSE))
          })

test_that("Argument 'combined_sample_name' creates error if vector with two 
          entries is provided.", {
  expect_error(combine_regions(data = test_data_filtered,
                               combined_sample_name = c("Consensus","Two"),
                               show_messages = FALSE))
  })

test_that("Argument 'combined_sample_name' creates error if 'NA' is 
          provided.", {
  expect_error(combine_regions(data = test_data_filtered,
                               combined_sample_name = NA,
                               show_messages = FALSE))
  })

### -----------------------------------------------------------------------###

test_that("Argument 'show_messages' creates no error if TRUE or FALSE
          value is provided.", {
            expect_no_error(combine_regions(data = test_data_filtered,
                                            show_messages = FALSE))
            expect_no_error(combine_regions(data = test_data_filtered,
                                            show_messages = TRUE))
          })

test_that("Argument 'show_messages' creates no error if non accepted
          value is provided.", {
            expect_error(combine_regions(data = test_data_filtered,
                                            show_messages = FaLSE))
          })

test_that("Argument 'show_messages' creates no error if non accepted
          value 'NA' is provided.", {
            expect_error(combine_regions(data = test_data_filtered,
                                            show_messages = NA))
          })


### -----------------------------------------------------------------------###
### Test input
### -----------------------------------------------------------------------###

test_that("Input data frame has the expected structure", {
  data <- test_data_filtered

  expect_equal(length(names(data)), 8)
  expect_identical(names(data), input_colnames)
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

### -----------------------------------------------------------------------###
### Test Output
### -----------------------------------------------------------------------###

test_that("Output data has the correct classes and structure", {
  expect_no_error(check_data_structure(test_data_combined))
})

test_that("Output data frame has correct colnames", {
  expect_true(any(colnames(data) %in% output_colnames))
})

test_that("Output data frame has correct class", {
  expect_identical(class(data)[2], "tbl")
})

test_that("Output data frame is expected values", {
  expect_identical(data$center[1], 500)
  expect_identical(sum(data$score), 1140)
})

### -----------------------------------------------------------------------###

test_that("Output data results has correct summit for 'nearest' peak", {
  data <- combine_regions(
    data = test_data_filtered,
    found_in_samples = 2,
    combined_center = "nearest",
    annotate_with_input_names = FALSE,
    combined_sample_name = "consensus_peak",
    show_messages = FALSE
  )
  
  expect_identical(data$center[7], 500)
  expect_identical(data$name[7], "consensus_peak|7")
})  

test_that("Output data results has correct summit for 'strongst' peak", {
  data <- combine_regions(
    data = test_data_filtered,
    found_in_samples = 2,
    combined_center = "strongest",
    annotate_with_input_names = FALSE,
    combined_sample_name = "consensus_peak",
    show_messages = FALSE
  )
  
  expect_identical(data$center[7], 600)
  expect_identical(data$name[7], "consensus_peak|7")
})

test_that("Output data results has correct summit for 'middle' peak", {
  
  data <- combine_regions(
    data =  test_data_filtered,
    found_in_samples = 2,
    combined_center = "middle",
    annotate_with_input_names = FALSE,
    combined_sample_name = "consensus_peak",
    show_messages = FALSE
  )
  
  expect_identical(data$center[7], 550)
  expect_identical(data$name[7], "consensus_peak|7")
  
})
### -----------------------------------------------------------------------###
