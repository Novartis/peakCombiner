#' Apply user-defined filtering options to genomic regions.
#'
#' @description
#' [peakCombiner::filterRegions] is an optional step that allows
#' inclusion or exclusion of genomic regions based on 4 different criteria:
#'
#' * Include regions by their chromosome names (optional).
#' * Exclude blacklisted regions (optional).
#' * Include regions above a given score (optional).
#' * Include top n regions per sample, ranked from highest to lowest score
#'   (optional).
#'
#' The accepted input is the PeakCombiner data frame is created from the
#' function [peakCombiner::prepareInputRegions].
#' Please see [peakCombiner::prepareInputRegions] for more details.
#'
#' The [peakCombiner::filterRegions] can be used multiple times on the same
#' data set, which allows a user to step-wise optimize selection criteria of
#' regions of interest.
#'
#' @details
#' This is an optional step which enables commonly-needed filtering steps to
#' focus in on the key genomic regions of interest. This can be useful
#' when there are many genomic regions identified in your peak-caller or
#' input BED files.
#'
#' [peakCombiner::filterRegions] can be used multiple times on the same data
#' set, allowing a user to select regions of interest using a step-wise
#' optimization approach.
#'
#' * `includeByChromosomeName` -   Retains only chromosomes that are in the
#'                                    provided vector. By not including
#'                                    mitochondrial, sex, or non-classical
#'                                    chromosomes, genomic regions found on
#'                                    these chromosomes can be removed. If set
#'                                    to 'NULL' (default), this step will be
#'                                    skipped (optional).
#' * `excludeByBlacklist` -         A GenomicRanges file, dataframe or tibble 
#'                                    can be provided listing the genomic 
#'                                    regions to remove (having `chrom` (
#'                                    `seqnames` for GenomicRanges) , `start`, 
#'                                    and `end` column names). If set to 'NULL' 
#'                                    (default), this step will be skipped 
#'                                    (optional).
#'                                    Please note that if there are not matching
#'                                    entries in the 'chrom' columns of input
#'                                    and blacklist, an information message is
#'                                    displayed. This can happend and does not
#'                                    cause any problems with the script.
#' * `includeAboveScoreCutoff` -   Single numeric value that defines the
#'                                    `score` threshold above which all genomic
#'                                    regions will be retained. The `score`
#'                                    column in the peakCombiner input data
#'                                    should be non-zero for this parameter to
#'                                    be used. It is populated by
#'                                    [peakCombiner::prepareInputRegions], and
#'                                    by default takes the value of -log10(FDR)
#'                                    if possible (e.g., using a .narrowPeak
#'                                    file from MACS2 as input). Importantly,
#'                                    applying this filter retains a variable
#'                                    number of genomic regions per sample, all
#'                                    having a score greater than the
#'                                    `includeAboveScoreCutoff` parameter. If
#'                                    set to 'NULL' (default), this step will
#'                                    be skipped (optional).
#' * `includeTopNScoring` -        Single numeric value that defines how many
#'                                    of the top scoring genomic regions (using
#'                                    the column `score`) are retained. All
#'                                    other genomic regions are discarded.
#'                                    Importantly, applying this filter retains
#'                                    `includeTopNScoring` regions per
#'                                    sample, which means that the minimum
#'                                    enrichment levels may vary between
#'                                    samples. Note that if multiple genomic
#'                                    regions have the same `score` cutoff
#'                                    value, then all of those genomic regions
#'                                    are included. In this case, the number of
#'                                    resulting regions retained may be a bit
#'                                    higher than the input parameter. If set to
#'                                    'NULL' (default), this step will be
#'                                    skipped (optional).
#'
#' @inheritParams centerExpandRegions
#'
#' @param includeByChromosomeName
#'          * 'NULL' (default) - No chromosome name filtering will be done.
#'          * Character vector that contains chromosomes names to be retained.
#'
#' @param excludeByBlacklist
#'          * 'NULL' (default) - No blacklist filtering will be done.
#'          * GenomicRanges object (default setup) or data frame/tibble with 
#'            columns `chrom`, `start`, and `end`.
#'
#' @param includeAboveScoreCutoff
#'          * 'NULL' (default) - No score filtering will be done.
#'          * Single numeric value that defines the `score` threshold above
#'            which all genomic regions will be retained. This results in
#'            variable number of sites per sample.
#'
#' @param includeTopNScoring
#'          * 'NULL' (default) - No score filtering will be done.
#'          * Single numeric value representing the number of genomic regions
#'            per sample to be retained. The genomic regions are selected from
#'            highest to lowest score, and if includeTopNScoring > number of
#'            regions, then no filtering is done.
#'            
#' @param outputFormat Character value to define format of output object. 
#'                      Accepted values are "GenomicRanges" (default), "tibble" 
#'                      or "data.frame".  
#'
#' @param showMessages Logical value of TRUE (default) or FALSE. Defines if
#'                      info messages are displayed or not.
#'
#' @return A tibble with the columns `chrom`, `start`, `end`, `name`, `score`,
#' `strand`, `center`, `sample_name`. The definitions of these columns are
#' described in full in the [peakCombiner::prepareInputRegions] Details.
#' Use as input for functions [peakCombiner::centerExpandRegions] and
#' [peakCombiner::combineRegions].
#'
#' @export
#'
#' @importFrom rlang .data
#' @import tidyr
#' @import here
#'
#' @examples
#'
#' # Load in and prepare a an accepted tibble
#' utils::data(syn_data_bed)
#'
#' data_prepared <- prepareInputRegions(
#'   data = syn_data_bed,
#'   outputFormat = "tibble",
#'   showMessages = TRUE
#' )
#'
#' # Here use options for all four filtering methods.
#'
#' filterRegions(
#'   data = data_prepared,
#'   includeByChromosomeName = c("chr1", "chr2", "chr4"),
#'   excludeByBlacklist = NULL,
#'   includeAboveScoreCutoff = 10,
#'   includeTopNScoring = 100,
#'   outputFormat = "tibble",
#'   showMessages = TRUE
#' )
#'
filterRegions <- function(data,
                           includeByChromosomeName = NULL,
                           excludeByBlacklist = NULL,
                           includeAboveScoreCutoff = NULL,
                           includeTopNScoring = NULL,
                           outputFormat = "GenomicRanges",
                           showMessages = TRUE) {
  
  
  ### -----------------------------------------------------------------------###
  ### Pre-Check up
  ### -----------------------------------------------------------------------###
  ## Check the validity of the peakCombiner input data format
  data <- checkDataStructure(
    data = data, 
    showMessages = showMessages
  )
  
  ### -----------------------------------------------------------------------###
  ### Show or hide messages
  ### -----------------------------------------------------------------------###
  
  if (!is.logical(showMessages)) {
    # show error message independent of parameter showMessages
    options("rlib_message_verbosity" = "default")
    
    cli::cli_abort(c(
      "x" = "Argument {.arg showMessages} has to be {.cls logical}."
    ))
  } else if (isTRUE(showMessages)) {
    options("rlib_message_verbosity" = "default")
  } else if (isFALSE(showMessages)) {
    options("rlib_message_verbosity" = "quiet")
  } else {
    # show error message independent of parameter showMessages
    options("rlib_message_verbosity" = "default")
    
    cli::cli_abort(c(
      "x" = "Argument {.arg showMessages} is a non-accepted {.cls logical}
      value.",
      "i" = "Argument {.arg showMessages} is {.val {showMessages}}."
    ))
  }
  
  ### -----------------------------------------------------------------------###
  ### Define parameters
  ### -----------------------------------------------------------------------###
  ##
  
  ### -----------------------------------------------------------------------###
  ### Check if GenomicRanges object contains only one genome
  ### -----------------------------------------------------------------------###
  
  if (inherits(data, "GRanges")) {
    cli::cli_inform(c(
      "i" = "Input data {.arg data} is a class {.cls GRanges}."
    ))    
    
    #' Extract info about genome from GenomicRanges file
    input_file_genome <- GenomeInfoDb::genome(data) |> unique()
    
    if (length(input_file_genome) > 1) {
      cli::cli_abort(c(
        "i" = "Input data {.arg data} is a class {.cls GRanges}.",
        "x" = "Input data {.arg data} has multiple assigned genomes.
        Input data has to have be from the same genome.",
        "i" = "Values of assigned genomes are: {.val {input_file_genome}}."
      ))
    }
    cli::cli_inform(c(
      "i" = "Input data {.arg data} assigned genomes is 
      {.val {input_file_genome}}."
    ))  
    
  }
  
  ##
  ## Pass data into new variable
  ### -----------------------------------------------------------------------###
  ### Figure out what kind of input data was entered by the user and
  ### load the initial data for follow-up quality checks
  ### -----------------------------------------------------------------------###
  
  required_colnames <- c(
    "chrom", "start", "end", "sample_name"
  )
  
  if (inherits(data, "GRanges")) {
    cli::cli_inform(c(
      "!" = "Provided input {.arg data} is a class {.cls GRanges} and will be
      converted to class {.cls tibble}.",
      ">" = "Start converting and preparing data."
    ))
    
    input_seqinfo <- GenomeInfoDb::seqinfo(data)
    
    data_filtered <-
      tibble::as_tibble(data) |>
      dplyr::rename(chrom = .data$seqnames) |>
      dplyr::mutate(
        start = as.numeric(.data$start),
        end = as.numeric(.data$end),
        strand = as.character(.data$strand)
      ) |>
      dplyr::mutate(strand = ifelse(.data$strand == "*", ".", .data$strand))
  } else if (all(required_colnames %in% colnames(data))) {
    cli::cli_inform(c(
      "i" = "Provide input {.arg data} is a {.cls data.frame} with three or four
      columns and paths to existing files.",
      ">" = "Start loading and preparing data."
    ))
    
    data_filtered <- data
    
  } else if (all(required_colnames %in% colnames(data))) {
    data_filtered <- data
    
    cli::cli_inform(c(
      "i" = "Provide input {.arg data} is a pre-loaded {.cls data.frame}  with
      the required column names.",
      ">" = "Start preparing data."
    ))
  } else {
    cli::cli_abort(c(
      "x" = "Provide input {.arg data} does not have the required format.",
      "!" = "Please check your column names in {.arg data}."
    ))
  }
  
  ### -----------------------------------------------------------------------###
  ### Check if output format is valid
  ### -----------------------------------------------------------------------###
  
  if (outputFormat %in% c("GenomicRanges", 
                           "GRanges", 
                           "tibble", 
                           "data.frame", 
                           "data.table")) {
    cli::cli_inform(c(
      "i" = "Argument {.arg outputFormat} is set to {.val {outputFormat}}."
    ))
  } else {
    cli::cli_abort(c(
      "x" = "Argument {.arg outputFormat} has to be one of the following
      values: {.val GenomicRanges}, {.val tibble}, or {.val data.frame}.",
      "i" = "Provided value is {.val {outputFormat}}."
    ))
  }
  
  ### -----------------------------------------------------------------------###
  ### Filter by chromosomes names
  ### -----------------------------------------------------------------------###

  data_filtered <-
    filterByChromosomeNames(
      data = data_filtered,
      includeByChromosomeName = includeByChromosomeName,
      showMessages = showMessages
    )

  ### -----------------------------------------------------------------------###
  ### Filter by blacklist
  ### -----------------------------------------------------------------------###

  data_filtered <-
    filterByBlacklist(
      data = data_filtered,
      excludeByBlacklist = excludeByBlacklist,
      showMessages = showMessages
    )

  ### -----------------------------------------------------------------------###
  ### Filter for significance values (score)
  ### -----------------------------------------------------------------------###

  data_filtered <-
    filterBySignificance(
      data = data_filtered,
      includeAboveScoreCutoff = includeAboveScoreCutoff,
      showMessages = showMessages
    )

  ### -----------------------------------------------------------------------###
  ### Select top peaks
  ### -----------------------------------------------------------------------###

  data_filtered <-
    filterByTopEnriched(
      data = data_filtered,
      includeTopNScoring = includeTopNScoring,
      showMessages = showMessages
    )

  ### -----------------------------------------------------------------------###
  ### Return data
  ### -----------------------------------------------------------------------###

  data_filtered <- data_filtered |>
    dplyr::relocate("strand", .after = "score") |>
    dplyr::mutate(strand = ifelse(.data$strand == "*", ".", .data$strand)) |>
    dplyr::ungroup()

  cli::cli_inform(c(
    "v" = "Filtered dataset will be returned."
  ))
  
  ### -----------------------------------------------------------------------###
  ### Adjust output format
  ### -----------------------------------------------------------------------###
  
  if (outputFormat %in% c("GenomicRanges", "GRanges")) {
    if(exists("input_seqinfo")) {
      cli::cli_inform(c(
        "i" = "Output format is set to {.val {outputFormat}}.",
        "i" = "Assigning input genome annotation to ouutput. ")
      )
      
      data_filtered <- 
        data_filtered |>
        GenomicRanges::makeGRangesFromDataFrame(
          keep.extra.columns = TRUE,
          seqinfo = input_seqinfo
        )
    } else{
      cli::cli_inform(c(
        "i" = "Output format is set to {.val {outputFormat}}.",
        "i" = "No input genome annotation assigned to ouutput. ")
      )
      data_filtered <- 
        data_filtered |>
        GenomicRanges::makeGRangesFromDataFrame(
          keep.extra.columns = TRUE
        )
    }
  } else if (outputFormat %in% c("tibble", "data.frame", "data.table")) {
    cli::cli_inform(c(
      "i" = "Output format is set to {.val tibble}."
    ))
  } else {
    cli::cli_abort(c(
      "x" = "Argument {.arg outputFormat} has to be one of the following
      values: {.val GenomicRanges}, {.val tibble}, or {.val data.frame}.",
      "i" = "Provided value is {.val {outputFormat}}."
    ))
  } 
  
  ### -----------------------------------------------------------------------###
  ### Set message display back to default
  ### -----------------------------------------------------------------------###

  if (isFALSE(showMessages)) {
    options("rlib_message_verbosity" = "default")
  }

  return(data_filtered)
}

