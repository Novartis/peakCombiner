#' Combine overlapping genomic regions from different samples to create
#' a single set of consensus genomic regions
#'
#' @description
#' [peakCombiner::combineRegions] is the main function of this package and
#' combines overlapping genomic regions from different samples to create
#' a single set of consensus genomic regions.
#'
#' The accepted input is the PeakCombiner data frame is created from the
#' function [peakCombiner::prepareInputRegions] and has optionally
#' already been centered and expanded and / or filtered using
#' [peakCombiner::centerExpandRegions] and [peakCombiner::filterRegions],
#' respectively.
#' Please see [peakCombiner::prepareInputRegions] for more details.
#'
#' @details
#' [peakCombiner::combineRegions] creates a set of consensus genomic regions by
#' combining overlapping genomic regions from different samples.
#' The general steps within this function are:
#'
#' * Identify overlapping genomic regions from the input samples
#' * Retain overlapping genomic regions that are found in at least
#'   `foundInSamples` samples. In this way, you can remove rare or
#'   sample-specific regions
#' * Note that overlapping genomic regions must contain at least one 'center'
#'   from its input sample regions to be considered a valid genomic region.
#' * As you can use the output data from this step again (e.g., to
#'   center and expand the new set of consensus regions), we must define
#'   the 'center', 'score', 'sample_name', and 'name' values for the new
#'   genomic regions. We do this as follows:
#'    + 'center' is defined by the `combinedCenter` parameter, which has three
#'       options.
#'          * `middle`        - the mathematical center of the new region
#'          * `strongest`     - the 'center' of the input region that has the
#'                              the highest 'score' of all overlapping input
#'                              regions
#'          * `nearest`       - the 'center' of the input region that is closest
#'                              to mean of the 'center's of all overlapping
#'                              input regions (default)
#'    + 'score' is the score of the genomic region from the sample whose
#'      'center's was used, or the mean of the 'score's if `middle` was selected
#'      for the `combinedCenter` parameter
#'    + 'sample_name' can be user defined (`combinedSampleName`) or is a
#'    concatenated string of all input 'sample_names' (default).
#'    + 'name' is created by combining 'sample_name' and row number to create a
#'    unique identifier for each newly created genomic region.
#'
#' Note, the output data.frame columns `sample_name`, `name` and `score`
#' will be updated.
#'
#' @param data        PeakCombiner data frame structure with required columns
#'                      named `chrom`, `start`, `end`, `name`,
#'                      `score`, `strand`, `center`, `sample_name`. Additional
#'                      columns will be dropped
#' @param foundInSamples  Only include genomic regions that are found
#'                            in at least `foundInSamples` **number**
#'                            of samples. If `foundInSamples` is a fraction
#'                            between 0 and 1, then only include genomic
#'                            regions that ar found in at least
#'                            `foundInSamples` **fraction** of samples.
#'                            Default value is 2.
#' @param combinedCenter   Defines how the column 'center' will be
#'                            populated for each genomic region in the output
#'                            data. Allowed options are
#'          * `middle`        - the mathematical center of the new region
#'          * `strongest`     - the 'center' of the input region that has the
#'                              the highest 'score' of all overlapping input
#'                              regions
#'          * `nearest`       - the 'center' of the input region that is closest
#'                              to mean of the 'center's of all overlapping
#'                              input regions (default)
#' @param annotateWithInputNames TRUE / FALSE (default). If TRUE, a new
#'                                    column named 'input_names' is created
#'                                    in the output data that is populated for
#'                                    each combined genomic region with the
#'                                    'name's of all contributing input regions.
#'                                    If the column 'input_names' already
#'                                    exists, it will be overwritten.
#' @param combinedSampleName Optionally defines how the column 'sample_name'
#'                               is populated for the output data.
#'                               If not used, then the default is to simply
#'                               concatenate all input
#'                               sample_names into a single comma-separated
#'                               string
#'            
#' @param outputFormat Character value to define format of output object. 
#'                      Accepted values are "GenomicRanges" (default), "tibble" 
#'                      or "data.frame".  
#'
#' @param showMessages Logical value of TRUE (default) or FALSE. Defines if
#'                      info messages are displayed or not.
#'
#' @return A tibble with the columns `chrom`, `start`, `end`, `name`, `score`,
#' `strand`, `center`, `sample_name`, and optionally `input_names`.
#' The definitions of these columns are
#' described in full in the Details below. Use as input for functions
#' [peakCombiner::centerExpandRegions] and [peakCombiner::filterRegions].
#'
#' @export
#'
#' @importFrom rlang .data
#' @import stringr
#' @import tidyr
#' @import here
#'
#' @examples
#' # Load in and prepare a an accepted tibble
#' utils::data(syn_data_bed)
#'
#' data_prepared <- prepareInputRegions(
#'   data = syn_data_bed,
#'   outputFormat = "tibble",
#'   showMessages = FALSE
#' )
#'
#' # Lets combine the input data by defining all potential option
#' combineRegions(
#'   data = data_prepared,
#'   foundInSamples = 2,
#'   combinedCenter = "nearest",
#'   annotateWithInputNames = TRUE,
#'   combinedSampleName = "consensus",
#'   outputFormat = "tibble",
#'   showMessages = TRUE
#' )
#'
combineRegions <- function(data,
                            foundInSamples = 2,
                            combinedCenter = "nearest",
                            annotateWithInputNames = FALSE,
                            combinedSampleName = NULL,
                            outputFormat = "GenomicRanges",
                            showMessages = TRUE) {
  ### -----------------------------------------------------------------------###
  ### Correct parameters & load needed variables
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
    # show error message independent of parameter showMessages
    options("rlib_message_verbosity" = "default")
    
    cli::cli_abort(c(
      "x" = "Argument {.arg outputFormat} has to be one of the following
      values: {.val GenomicRanges}, {.val tibble}, or {.val data.frame}.",
      "i" = "Provided value is {.val {outputFormat}}."
    ))
  }
  
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
    # show error independend of showMessages
    options("rlib_message_verbosity" = "default")
    
    cli::cli_abort(c(
      "x" = "Provide input {.arg data} does not have the required format.",
      "!" = "Please check your column names in {.arg data}."
    ))
  }
  
  ##
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
  ### Pre-Check up
  ### -----------------------------------------------------------------------###
  ## Check the validity of the peakCombiner input data format

  data <- checkDataStructure(
    data = data
  )

  ### -----------------------------------------------------------------------###
  ### Combine peaks - Disjoin & Filter
  ### -----------------------------------------------------------------------###
  ## 1: Do a disjoin to separate the peaks and filter based on foundInSamples
  data_disjoin <- crDisjoinFilter(
    data = data,
    foundInSamples = foundInSamples
  )

  ### -----------------------------------------------------------------------###
  ### Combine peaks - Reduce
  ### -----------------------------------------------------------------------###
  ## 2: Reduce the disjoined data and prepare combined table
  data_reduce <- crReduce(
    data = data_disjoin
  )

  ### -----------------------------------------------------------------------###
  ### Combine peaks - Overlap with summit
  ### -----------------------------------------------------------------------###
  ## 3: Remove false positive peaks without summit
  data_overlap_summit <- crOverlapWithSummits(
    data = data_reduce,
    input = data
  )

  ### -----------------------------------------------------------------------###
  ### Combine peaks - Link to best summit
  ### -----------------------------------------------------------------------###
  ## 4: Identify top enriched summit for new defined peaks
  data_combined_with_summit <- crAddSummit(
    data = data_overlap_summit,
    input = data,
    combinedCenter = combinedCenter,
    annotateWithInputNames = annotateWithInputNames,
    combinedSampleName = combinedSampleName
  )


  data_combined_with_summit <- data_combined_with_summit |>
    dplyr::relocate("strand", .after = "score") |>
    dplyr::mutate(strand = ifelse(.data$strand == "*", ".", .data$strand)) |>
    dplyr::ungroup()


  ### -----------------------------------------------------------------------###
  ### Combine peaks - Return data frame
  ### -----------------------------------------------------------------------###

  cli::cli_inform(c(
    "v" = "Genomic regions were successfully combined.",
    " " = " "
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
      
      data_combined_with_summit <- 
        data_combined_with_summit |>
        GenomicRanges::makeGRangesFromDataFrame(
          keep.extra.columns = TRUE,
          seqinfo = input_seqinfo
        )
    } else{
      cli::cli_inform(c(
        "i" = "Output format is set to {.val {outputFormat}}.",
        "i" = "No input genome annotation assigned to ouutput. ")
      )
      data_combined_with_summit <- 
        data_combined_with_summit |>
        GenomicRanges::makeGRangesFromDataFrame(
          keep.extra.columns = TRUE
        )
    }
  } else if (outputFormat %in% c("tibble", "data.frame", "data.table")) {
    cli::cli_inform(c(
      "i" = "Output format is set to {.val tibble}."
    ))
  } else {
    # show error message independent of parameter showMessages
    options("rlib_message_verbosity" = "default")
    
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

  return(data_combined_with_summit)
}
