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
test_data_prepared <- peakCombiner::prepareInputRegions(
  data = test_data,
  outputFormat = "tibble"
)
##
test_data_center_expand <- peakCombiner::centerExpandRegions(
  data = test_data_prepared,
  centerBy = "center_column",
  outputFormat = "tibble",
  expandBy = 9e0,
  showMessages = FALSE
)

test_data_center_expand <- peakCombiner::centerExpandRegions(
  data = test_data_prepared,
  centerBy = "midpoint",
  outputFormat = "tibble",
  genome = "hg38",
  expandBy = 90,
  showMessages = F
)


restult_colnames <- colnames(test_data_center_expand)
##
test_data_filtered <- peakCombiner::filterRegions(
  data = test_data_center_expand,
  excludeByBlacklist = NULL,
  includeByChromosomeName = NULL,
  includeAboveScoreCutoff = NULL,
  includeTopNScoring = NULL,
  outputFormat = "tibble",
  showMessages = FALSE
)
##
test_data_combined <- peakCombiner::combineRegions(
  data = test_data_filtered,
  foundInSamples = 2,
  combinedCenter = "nearest",
  annotateWithInputNames = TRUE,
  combinedSampleName = "combined",
  outputFormat = "tibble",
  showMessages = FALSE
)
##
test_data_combined_ce <- peakCombiner::centerExpandRegions(
  data = test_data_combined,
  centerBy = "center_column",
  outputFormat = "tibble",
  expandBy = NULL,
  showMessages = FALSE
)


### -----------------------------------------------------------------------###
### Test genome
### -----------------------------------------------------------------------###

testthat::test_that("Test if function works with genome input", {
  testthat::expect_no_error(peakCombiner::centerExpandRegions(
    data = test_data_prepared,
    centerBy = "center_column",
    expandBy = NULL,
    outputFormat = "tibble",
    genome = "hg38",
    trim_start = FALSE,
    showMessages = FALSE
  ))
})

testthat::test_that("Test if function works with genome input & trimming", {
  testthat::expect_no_error(peakCombiner::centerExpandRegions(
    data = test_data_prepared,
    centerBy = "center_column",
    expandBy = NULL,
    outputFormat = "tibble",
    genome = "hg38",
    trim_start = TRUE,
    showMessages = FALSE
  ))
})

testthat::test_that("Test if function fails with wrong genome", {
  testthat::expect_error(peakCombiner::centerExpandRegions(
    data = test_data_prepared,
    centerBy = "center_column",
    expandBy = NULL,
    outputFormat = "tibble",
    genome = "HG39",
    trim_start = TRUE,
    showMessages = FALSE
  ))
})

testthat::test_that("Test if function fails with wrong genome", {
  testthat::expect_error(peakCombiner::centerExpandRegions(
    data = test_data_prepared,
    centerBy = "center_column",
    expandBy = NULL,
    outputFormat = "tibble",
    genome = "mm38",
    trim_start = TRUE,
    showMessages = FALSE
  ))
})

testthat::test_that("Test if function fails with wrong trim", {
  testthat::expect_error(peakCombiner::centerExpandRegions(
    data = test_data_prepared,
    centerBy = "center_column",
    expandBy = NULL,
    outputFormat = "tibble",
    genome = "hg38",
    trim_start = 10,
    showMessages = FALSE
  ))
})

testthat::test_that("Test if function fails with wrong trim", {
  testthat::expect_no_error(peakCombiner::centerExpandRegions(
    data = test_data_prepared,
    centerBy = "center_column",
    expandBy = NULL,
    outputFormat = "tibble",
    genome = NA,
    trim_start = TRUE,
    showMessages = FALSE
  ))
})


### -----------------------------------------------------------------------###
### Test input
### -----------------------------------------------------------------------###

testthat::test_that("Test if function works with pre-combined input", {
  testthat::expect_no_error(peakCombiner::centerExpandRegions(
    data = test_data_prepared,
    centerBy = "center_column",
    expandBy = NULL,
    genome = NA,
    trim_start = FALSE,
    outputFormat = "tibble"
  ))
})

testthat::test_that("Test if function works with post-combined input", {
  testthat::expect_no_error(peakCombiner::centerExpandRegions(
    data = test_data_combined,
    centerBy = "center_column",
    expandBy = NULL,
    genome = NA,
    trim_start = FALSE,
    outputFormat = "tibble"
  ))
})

### -----------------------------------------------------------------------###

testthat::test_that("Required input data has the expected structure", {
  data <- test_data_prepared
  
  testthat::expect_equal(length(names(data)), 8)
  testthat::expect_identical(names(data), input_colnames_pre)
  testthat::expect_true(is.character(data$chrom))
  testthat::expect_true(is.numeric(data$start))
  testthat::expect_true(is.numeric(data$end))
  testthat::expect_true(is.character(data$name))
  testthat::expect_true(is.numeric(data$score))
  testthat::expect_true(is.character(data$strand))
  testthat::expect_true(is.numeric(data$center))
  testthat::expect_true(is.character(data$sample_name))
  testthat::expect_true(sum(stringr::str_detect(data$name, "|")) > 0)
})

testthat::test_that("Required input data has the expected structure", {
  data <- test_data_combined
  
  testthat::expect_equal(length(names(data)), 9)
  testthat::expect_identical(names(data), input_colnames_post)
  testthat::expect_true(is.character(data$chrom))
  testthat::expect_true(is.numeric(data$start))
  testthat::expect_true(is.numeric(data$end))
  testthat::expect_true(is.character(data$name))
  testthat::expect_true(is.numeric(data$score))
  testthat::expect_true(is.character(data$strand))
  testthat::expect_true(is.numeric(data$center))
  testthat::expect_true(sum(stringr::str_detect(data$input_names, "|")) > 0)
})

### -----------------------------------------------------------------------###

testthat:: test_that("Required paramter 'centerBy' has the expected structure/value", {
  testthat::expect_no_error(peakCombiner::centerExpandRegions(
    data = test_data_prepared,
    centerBy = "center_Column",
    expandBy = NULL,
    genome = NA,
    trim_start = FALSE,
    outputFormat = "tibble"
  ))
  testthat::expect_error(peakCombiner::centerExpandRegions(
    data = test_data_prepared,
    centerBy = c("center_column", "calculated_value"),
    expandBy = NULL,
    genome = NA,
    trim_start = FALSE,
    outputFormat = "tibble"
  ), "centerBy")
  testthat::expect_error(peakCombiner::centerExpandRegions(
    data = test_data_prepared,
    centerBy = "nonexisting",
    expandBy = NULL,
    genome = NA,
    trim_start = FALSE,
    outputFormat = "tibble"
  ), "centerBy")
  testthat::expect_error(peakCombiner::centerExpandRegions(
    data = test_data_prepared,
    centerBy = NULL,
    expandBy = NULL,
    genome = NA,
    trim_start = FALSE,
    outputFormat = "tibble"
  ), "centerBy")
  testthat::expect_error(peakCombiner::centerExpandRegions(
    data = test_data_prepared,
    centerBy = NA,
    expandBy = NULL,
    genome = NA,
    trim_start = FALSE,
    outputFormat = "tibble"
  ), "centerBy")
})

### -----------------------------------------------------------------------###

testthat::test_that("Required paramter expandBy has the expected structure/value", {
  testthat::expect_no_error(peakCombiner::centerExpandRegions(
    data = test_data_prepared,
    centerBy = "center_column",
    expandBy = NULL,
    genome = NA,
    trim_start = FALSE,
    outputFormat = "tibble"
  ))
  testthat::expect_error(peakCombiner::centerExpandRegions(
    data = test_data_prepared,
    centerBy = "column_value",
    expandBy = NA,
    genome = NA,
    trim_start = FALSE,
    outputFormat = "tibble"
  ), )
  testthat::expect_error(peakCombiner::centerExpandRegions(
    data = test_data_prepared,
    centerBy = "column_value",
    expandBy = c(1, 2, 3),
    genome = NA,
    trim_start = FALSE,
    outputFormat = "tibble"
  ), )
  testthat::expect_error(peakCombiner::centerExpandRegions(
    data = test_data_prepared,
    centerBy = "column_value",
    expandBy = "nonexisting",
    genome = NA,
    trim_start = FALSE,
    outputFormat = "tibble"
  ), )
})

### -----------------------------------------------------------------------###
### Test Output
### -----------------------------------------------------------------------###

testthat:: test_that("Output data frame is correct for pre-combined", {
  data <- test_data_center_expand
  
  testthat::expect_setequal(colnames(data), output_colnames_pre)
  testthat::expect_equal(ncol(data), 8)
  
  testthat::expect_identical(class(data)[2], "tbl")
  
  testthat::expect_true(is.character(data$chrom))
  testthat::expect_true(is.numeric(data$start))
  testthat::expect_true(is.numeric(data$end))
  testthat::expect_true(is.character(data$name))
  testthat::expect_true(is.numeric(data$score))
  testthat::expect_true(is.character(data$strand))
  testthat::expect_true(is.numeric(data$center))
  testthat::expect_true(is.character(data$sample_name))
  
  testthat::expect_equal(mean(data$center), 2954.761905)
  testthat::expect_identical(nrow(data), as.integer(42))
  testthat::expect_identical(data$start[1], 559L)
})

test_that("Output data frame is correct for post-combined", {
  data <- test_data_combined
  
  testthat::expect_setequal(colnames(data), output_colnames_post)
  testthat::expect_equal(ncol(data), 9)
  
  testthat::expect_identical(class(data)[2], "tbl")
  
  testthat::expect_true(is.character(data$chrom))
  testthat::expect_true(is.numeric(data$start))
  testthat::expect_true(is.numeric(data$end))
  testthat::expect_true(is.character(data$name))
  testthat::expect_true(is.numeric(data$score))
  testthat::expect_true(is.character(data$strand))
  testthat::expect_true(is.numeric(data$center))
  testthat::expect_equal(mean(data$center), 3318.750)
  testthat::expect_identical(nrow(data), as.integer(8))
  testthat::expect_identical(data$start[1], 459)
  testthat::expect_identical(data$end[1], 738)
})

testthat:: test_that("Output data frame is correct for data_prepared", {
  ##
  data <- test_data_prepared
  result <- peakCombiner::centerExpandRegions(
    data = data,
    centerBy = "center_column",
    expandBy = NULL,
    genome = NA,
    trim_start = FALSE,
    outputFormat = "tibble"
  )
  ##
  testthat::expect_no_error(peakCombiner::centerExpandRegions(
    data = data,
    centerBy = "center_column",
    expandBy = NULL,
    genome = NA,
    trim_start = FALSE,
    outputFormat = "tibble"
  ))
  ##
  testthat::expect_identical(nrow(result), 51L)
})

##
testthat:: test_that("Output data frame is correct for data_center_expand", {
  ##
  data <- test_data_center_expand
  result <- peakCombiner::centerExpandRegions(
    data = data,
    centerBy = "center_column",
    expandBy = NULL,
    genome = NA,
    trim_start = FALSE,
    outputFormat = "tibble"
  )
  ##
  testthat::expect_no_error(peakCombiner::centerExpandRegions(
    data = data,
    centerBy = "center_column",
    expandBy = NULL,
    genome = NA,
    trim_start = FALSE,
    outputFormat = "tibble"
  ))
  ##
  testthat::expect_identical(nrow(result), 42L)
})
##
testthat:: test_that("Output data frame is correct for data_filtered", {
  ##
  data <- test_data_filtered
  result <- peakCombiner::centerExpandRegions(
    data = data,
    centerBy = "center_column",
    expandBy = NULL,
    genome = NA,
    trim_start = FALSE,
    outputFormat = "tibble"
  )
  ##
  testthat::expect_no_error(peakCombiner::centerExpandRegions(
    data = data,
    centerBy = "center_column",
    expandBy = NULL,
    genome = NA,
    trim_start = FALSE,
    outputFormat = "tibble"
  ))
  ##
  testthat::expect_identical(nrow(result), 42L)
})
##
testthat:: test_that("Output data frame is correct for data_combined", {
  ##
  data <- test_data_combined
  result <- peakCombiner::centerExpandRegions(
    data = data,
    centerBy = "midpoint",
    expandBy = NULL,
    genome = NA,
    trim_start = FALSE,
    outputFormat = "tibble",
    showMessages = FALSE
  )
  testthat::expect_identical(nrow(result), 8L)
  
  ##
  testthat::expect_no_error(peakCombiner::centerExpandRegions(
    data = data,
    centerBy = "midpoint",
    expandBy = NULL,
    genome = NA,
    trim_start = FALSE,
    outputFormat = "tibble",
    showMessages = FALSE
  ))
  ##
})
##
### -----------------------------------------------------------------------###
##
testthat:: test_that("Check if default output is GRanges", {
  testthat::expect_no_error(peakCombiner::centerExpandRegions(
      data = test_data_prepared,
      genome = "hg38",
      showMessages = FALSE) 
    |> inherits("GRanges")
    )
})

#testthat:: test_that("Check if genome is hg38 in output", {
#  testthat::expect_no_error(
#    testthat::expect_identical(Seqinfo::genome(
#    peakCombiner::centerExpandRegions(
#    data = test_data_prepared,
#    genome = "hg38",
#    showMessages = FALSE) 
#    ) |> unique(), "hg38")
#  )
#})

testthat:: test_that("Using a not supported genome", {
  testthat::expect_error(
   peakCombiner::centerExpandRegions(
      data = test_data_prepared,
      genome = "mm9",
      showMessages = FALSE) 
    )
})

testthat:: test_that("Using a not wrong spelling of the genome", {
  testthat::expect_error(
    peakCombiner::centerExpandRegions(
      data = test_data_prepared,
      genome = "MM10",
      showMessages = FALSE) 
  )
})

testthat:: test_that("Using a not wrong spelling of the genome", {
  testthat::expect_error(
    peakCombiner::centerExpandRegions(
      data = test_data_prepared,
      genome = c("mm10", "hg38"),
      showMessages = FALSE) 
  )
})

