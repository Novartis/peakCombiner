#' Combine overlapping genomic regions from different samples to create
#' a single set of consensus genomic regions
#'
#' @description
#' [peakCombiner::combine_regions] is the main function of this package and
#' combines overlapping genomic regions from different samples to create
#' a single set of consensus genomic regions.
#'
#' The accepted input is the PeakCombiner data frame is created from the
#' function [peakCombiner::prepare_input_regions] and has optionally
#' already been centered and expanded and / or filtered using
#' [peakCombiner::center_expand_regions] and [peakCombiner::filter_regions],
#' respectively.
#' Please see [peakCombiner::prepare_input_regions] for more details.
#'
#' @details
#' [peakCombiner::combine_regions] creates a set of consensus genomic regions by
#' combining overlapping genomic regions from different samples.
#' The general steps within this function are:
#'
#' * Identify overlapping genomic regions from the input samples
#' * Retain overlapping genomic regions that are found in at least
#'   `found_in_samples` samples. In this way, you can remove rare or
#'   sample-specific regions
#' * Note that overlapping genomic regions must contain at least one 'center'
#'   from its input sample regions to be considered a valid genomic region.
#' * As you can use the output data from this step again (e.g., to
#'   center and expand the new set of consensus regions), we must define
#'   the 'center', 'score', 'sample_name', and 'name' values for the new
#'   genomic regions. We do this as follows:
#'    + 'center' is defined by the `combined_center` parameter, which has three
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
#'      for the `combined_center` parameter
#'    + 'sample_name' can be user defined (`combined_sample_name`) or is a
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
#' @param found_in_samples  Only include genomic regions that are found
#'                            in at least `found_in_samples` **number**
#'                            of samples. If `found_in_samples` is a fraction
#'                            between 0 and 1, then only include genomic
#'                            regions that ar found in at least
#'                            `found_in_samples` **fraction** of samples.
#'                            Default value is 2.
#' @param combined_center   Defines how the column 'center' will be
#'                            populated for each genomic region in the output
#'                            data. Allowed options are
#'          * `middle`        - the mathematical center of the new region
#'          * `strongest`     - the 'center' of the input region that has the
#'                              the highest 'score' of all overlapping input
#'                              regions
#'          * `nearest`       - the 'center' of the input region that is closest
#'                              to mean of the 'center's of all overlapping
#'                              input regions (default)
#' @param annotate_with_input_names TRUE / FALSE (default). If TRUE, a new
#'                                    column named 'input_names' is created
#'                                    in the output data that is populated for
#'                                    each combined genomic region with the
#'                                    'name's of all contributing input regions.
#'                                    If the column 'input_names' already
#'                                    exists, it will be overwritten.
#' @param combined_sample_name Optionally defines how the column 'sample_name'
#'                               is populated for the output data.
#'                               If not used, then the default is to simply
#'                               concatenate all input
#'                               sample_names into a single comma-separated
#'                               string
#'
#' @param show_messages Logical value of TRUE (default) or FALSE. Defines if
#'                      info messages are displayed or not.
#'
#' @return A tibble with the columns `chrom`, `start`, `end`, `name`, `score`,
#' `strand`, `center`, `sample_name`, and optionally `input_names`.
#' The definitions of these columns are
#' described in full in the Details below. Use as input for functions
#' [peakCombiner::center_expand_regions] and [peakCombiner::filter_regions].
#'
#' @export
#'
#' @importFrom rlang .data
#' @import tidyr
#' @import here
#'
#' @examples
#' # Load in and prepare a an accepted tibble
#' input_data <- peakCombiner::syn_data_bed
#' input_data
#'
#' data_prepared <- prepare_input_regions(
#'   data = input_data,
#'   show_messages = FALSE
#' )
#'
#' # Lets combine the input data by defining all potential option
#' combine_regions(
#'   data = data_prepared,
#'   found_in_samples = 2,
#'   combined_center = "nearest",
#'   annotate_with_input_names = TRUE,
#'   combined_sample_name = "consensus",
#'   show_messages = TRUE
#' )
#'
combine_regions <- function(data,
                            found_in_samples = 2,
                            combined_center = "nearest",
                            annotate_with_input_names = FALSE,
                            combined_sample_name = NULL,
                            show_messages = TRUE) {
  ### -----------------------------------------------------------------------###
  ### Correct parameters & load needed variables
  ### -----------------------------------------------------------------------###

  ### -----------------------------------------------------------------------###
  ### Show or hide messages
  ### -----------------------------------------------------------------------###

  if (!is.logical(show_messages)) {
    # show error message independent of parameter show_messages
    options("rlib_message_verbosity" = "default")

    cli::cli_abort(c(
      "x" = "Argument {.arg show_messages} has to be {.cls logical}."
    ))
  } else if (isTRUE(show_messages)) {
    options("rlib_message_verbosity" = "default")
  } else if (isFALSE(show_messages)) {
    options("rlib_message_verbosity" = "quiet")
  } else {
    # show error message independent of parameter show_messages
    options("rlib_message_verbosity" = "default")

    cli::cli_abort(c(
      "x" = "Argument {.arg show_messages} is a non-accepted {.cls logical}
      value.",
      "i" = "Argument {.arg show_messages} is {.val {show_messages}}."
    ))
  }

  ### -----------------------------------------------------------------------###
  ### Pre-Check up
  ### -----------------------------------------------------------------------###
  ## Check the validity of the peakCombiner input data format

  data <- check_data_structure(
    data = data
  )

  ### -----------------------------------------------------------------------###
  ### Combine peaks - Disjoin & Filter
  ### -----------------------------------------------------------------------###
  ## 1: Do a disjoin to separate the peaks and filter based on found_in_samples
  data_disjoin <- cr_disjoin_filter(
    data = data,
    found_in_samples = found_in_samples
  )

  ### -----------------------------------------------------------------------###
  ### Combine peaks - Reduce
  ### -----------------------------------------------------------------------###
  ## 2: Reduce the disjoined data and prepare combined table
  data_reduce <- cr_reduce(
    data = data_disjoin
  )

  ### -----------------------------------------------------------------------###
  ### Combine peaks - Overlap with summit
  ### -----------------------------------------------------------------------###
  ## 3: Remove false positive peaks without summit
  data_overlap_summit <- cr_overlap_with_summits(
    data = data_reduce,
    input = data
  )

  ### -----------------------------------------------------------------------###
  ### Combine peaks - Link to best summit
  ### -----------------------------------------------------------------------###
  ## 4: Identify top enriched summit for new defined peaks
  data_combined_with_summit <- cr_add_summit(
    data = data_overlap_summit,
    input = data,
    combined_center = combined_center,
    annotate_with_input_names = annotate_with_input_names,
    combined_sample_name = combined_sample_name
  )


  data_combined_with_summit <- data_combined_with_summit |>
    dplyr::relocate(.data$strand, .after = .data$score) |>
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
  ### Set message display back to default
  ### -----------------------------------------------------------------------###

  if (isFALSE(show_messages)) {
    options("rlib_message_verbosity" = "default")
  }

  return(data_combined_with_summit)
}
