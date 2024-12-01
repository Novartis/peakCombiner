#' Synthetic sample sheet to load example data with peakCombiner
#'
#' Synthetic example sample sheet as tibble with columns "sample_name",
#' "file_path", "file_format", and "score_colname".
#'
#'
#' @format `syn_sample_sheet` A tibble with 6 rows and 4 columns.
#'
#' @source Created for R package peakCombiner.
#' @usage data(syn_sample_sheet)
"syn_sample_sheet"

#' Synthetic file with blacklisted regions for peakCombiner
#'
#' Synthetic example blacklisted regions file as tibble with columns "chrom",
#' "start", and "end".
#'
#' @format `syn_blacklist` A tibble with 2 rows and 3 columns:
#'
#' @source Created for R package peakCombiner.
#' @usage data(syn_blacklist)
"syn_blacklist"

#' Synthetic data set of genomic coordinates and meta data columns as tibble
#'
#' Synthetic example data set as tibble with columns "chrom", "start", "end",
#' "name", "score", "strand" , "center", and "sample_name".
#'
#'
#' @format `syn_data_tibble` A tibble with 55 rows and 8 columns:
#'
#' @source Created for R package peakCombiner.
#' @usage data(syn_data_tibble)
"syn_data_tibble"

#' Synthetic data set of genomic coordinates and meta data columns
#'
#' Synthetic example data set from GenomicRanges object with columns "seqnames",
#'  "start", "end", "width", "strand", "score", "center", and "sample_name".
#'
#'
#' @format `syn_data_granges` A data frame with 55 rows and 8 columns:
#'
#' @source Created for R package peakCombiner.
#' @usage data(syn_data_granges)
"syn_data_granges"

#' Synthetic data set of genomic coordinates and meta data columns
#'
#' Synthetic example data set as minimal required input file with columns
#' "chrom",  "start", "end", and "sample_name".
#'
#'
#' @format `syn_data_bed` A tibble with 55 rows and 4 columns:
#'
#' @source Created for R package peakCombiner.
#' @usage data(syn_data_bed)
"syn_data_bed"

#' Synthetic data set of genomic coordinates and meta data columns filtered for
#' control rep 1 sample
#'
#' Synthetic example data set as minimal required input file with columns
#' "chrom",  "start", "end", "score", "strand", and "center".
#'
#'
#' @format `syn_data_control01` A tibble with 11 rows and 6 columns:
#'
#' @source Created for R package peakCombiner.
#' @usage data(syn_data_control01)
"syn_data_control01"
#' Synthetic data set of genomic coordinates and meta data columns filtered for
#' treatment rep 1 sample
#'
#' Synthetic example data set as minimal required input file with columns
#' "chrom",  "start", "end", "score", "strand", and "center".
#'
#'
#' @format `syn_data_treatment01` A tibble with 10 rows and 6 columns:
#'
#' @source Created for R package peakCombiner.
#' @usage data(syn_data_treatment01)
"syn_data_treatment01"

#' Synthetic data set for control rep 1 sample in narrowPeak file format
#'
#' Synthetic example data set as minimal required input file with columns
#' "chrom",  "start", "end", "name", "score", "strand", "signalValue",
#' "pValue", "qValue" and "peak".
#'
#'
#' @format `syn_control_rep1_narrowPeak` A tibble with 11 rows and 6 columns:
#'
#' @source Created for R package peakCombiner.
#' @usage data(syn_control_rep1_narrowPeak)
"syn_control_rep1_narrowPeak"

#' Synthetic data set for treatment rep 1 sample in narrowPeak file format
#'
#' Synthetic example data set as minimal required input file with columns
#' "chrom",  "start", "end", "name", "score", "strand", "signalValue",
#' "pValue", "qValue" and "peak".
#'
#'
#' @format `syn_treatment_rep1_narrowPeak` A tibble with 11 rows and 6 columns:
#'
#' @source Created for R package peakCombiner.
#' @usage data(syn_treatment_rep1_narrowPeak)
"syn_treatment_rep1_narrowPeak"

#' Blacklisted regions from ENCODE for human hg38
#'
#' BED file format with blacklisted regions for human annotation hg38 with
#' column named "chrom", "start", and "end".
#'
#'
#' @format `blacklist_hg38` A tibble with 910 rows and 3 columns:
#'
#' @source Downloaded from ENCODE https://www.encodeproject.org/files/ENCFF356LFX/
#' @usage data(blacklist_hg38)
"blacklist_hg38"

#' Blacklisted regions from ENCODE for mouse mm10
#'
#' BED file format with blacklisted regions for mouse annotation mm10 with
#' column named "chrom", "start", and "end".
#'
#'
#' @format `blacklist_mm10` A tibble with 164 rows and 3 columns:
#'
#' @source Downloaded from ENCODE https://www.encodeproject.org/files/ENCFF547MET/
#' @usage data(blacklist_mm10)
"blacklist_mm10"
