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
test_data_prepared <- peakCombiner::prepareInputRegions(
  data = test_data,
  outputFormat = "tibble",
  showMessages = FALSE
)
##
test_data_center_expand <- peakCombiner::centerExpandRegions(
  data = test_data_prepared,
  centerBy = "center_column",
  expandBy = NULL,
  outputFormat = "tibble"
)
##
test_data_filtered <- peakCombiner:::filterByChromosomeNames(
  data = test_data_center_expand,
  includeByChromosomeName = c("chr1", "chr10", "chr42")
)
##
input_colnames <- colnames(test_data_filtered)
##
blacklist <- backlist <- tibble::tibble(chrom = c("chr1"),
                                start = c(100),
                                end = c(1000))
backlist

##
test_data_filtered_bl <- peakCombiner:::filterByBlacklist(
  data = test_data_filtered,
  excludeByBlacklist = blacklist
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
  expect_no_error(peakCombiner:::filterByBlacklist(
    data = test_data_filtered,
    excludeByBlacklist = blacklist
  ))
  expect_no_error(peakCombiner:::filterByBlacklist(
    data = test_data_filtered,
    excludeByBlacklist = NULL
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
test_that("Required parameter 'filterByBlacklist' has expected structure", {
  expect_no_error(peakCombiner:::filterByBlacklist(
    data = test_data_filtered,
    excludeByBlacklist = NULL
  ))
  ##
  expect_error(peakCombiner:::filterByBlacklist(
    data = test_data_filtered,
    excludeByBlacklist = "HG38"
  ))
  expect_error(peakCombiner:::filterByBlacklist(
    data = test_data_filtered,
    excludeByBlacklist = blacklist[1:2]
  ))
})
##
test_that("Required parameter 'filterByBlacklist' is GRanges object", {
  blacklist_gr <- blacklist |>
    GenomicRanges::makeGRangesFromDataFrame(
      keep.extra.columns = TRUE, 
    ) 
  expect_no_error(peakCombiner:::filterByBlacklist(
    data = test_data_filtered,
    excludeByBlacklist = blacklist_gr
  ))
  expect_no_error(peakCombiner:::filterByBlacklist(
    data = test_data_filtered,
    excludeByBlacklist = blacklist
  ))
})



##
### -----------------------------------------------------------------------###
##
test_that("For 'filterByBlacklist' providing blacklist with different
          names", {
  blacklist2 <- blacklist
  colnames(blacklist2) <- c("CHROM", "start", "end")
  ##
  expect_no_error(peakCombiner:::filterByBlacklist(
    data = test_data_filtered,
    excludeByBlacklist = blacklist2
  ))
  ##
  colnames(blacklist2) <- c("seqnames", "start", "end")
  ##
  expect_error(peakCombiner:::filterByBlacklist(
    data = test_data_filtered,
    excludeByBlacklist = blacklist2
  ))
})
##
### -----------------------------------------------------------------------###
##
test_that("Wrong input for excludeByBlacklist for 'filterByBlacklist'", {
  expect_error(peakCombiner:::filterByBlacklist(
    data = test_data_filtered,
    excludeByBlacklist = "mm38"
  ))
  expect_error(peakCombiner:::filterByBlacklist(
    data = test_data_filtered,
    excludeByBlacklist = hg38
  ))
  expect_error(peakCombiner:::filterByBlacklist(
    data = test_data_filtered,
    excludeByBlacklist = 1
  ))
  expect_error(peakCombiner:::filterByBlacklist(
    data = test_data_filtered,
    excludeByBlacklist = c(1, 2)
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
  expect_equal(round(mean(data$center), 0), 3883.0)
  expect_identical(nrow(data), 30L)
  expect_identical(data$start[1], 3450L)
})
##
### -----------------------------------------------------------------------###
