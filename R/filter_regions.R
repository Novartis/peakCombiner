#' Apply user-defined filtering options to genomic regions.
#'
#' @description
#' [peakCombiner::filter_regions] is an optional step that allows
#' inclusion or exclusion of genomic regions based on 4 different criteria:
#'
#' * Include regions by their chromosome names (optional).
#' * Exclude blacklisted regions (optional).
#' * Include regions above a given score (optional).
#' * Include top n regions per sample, ranked from highest to lowest score
#'   (optional).
#'
#' The accepted input is the PeakCombiner data frame is created from the
#' function [peakCombiner::prepare_input_regions].
#' Please see [peakCombiner::prepare_input_regions] for more details.
#'
#' The [peakCombiner::filter_regions] can be used multiple times on the same
#' data set, which allows a user to step-wise optimize selection criteria of
#' regions of interest.
#'
#' @details
#' This is an optional step which enables commonly-needed filtering steps to
#' focus in on the key genomic regions of interest. This can be useful
#' when there are many genomic regions identified in your peak-caller or
#' input BED files.
#'
#' [peakCombiner::filter_regions] can be used multiple times on the same data
#' set, allowing a user to select regions of interest using a step-wise
#' optimization approach.
#'
#' * `include_by_chromosome_name` -   Retains only chromosomes that are in the
#'                                    provided vector. By not including
#'                                    mitochondrial, sex, or non-classical
#'                                    chromosomes, genomic regions found on
#'                                    these chromosomes can be removed. If set
#'                                    to 'NULL' (default), this step will be
#'                                    skipped (optional).
#' * `exclude_by_blacklist` -         Remove ENCODE-annotated blacklisted
#'                                    regions for either human (
#'                      [hg38](https://www.encodeproject.org/files/ENCFF356LFX/)
#'                                    only) or mouse (
#'                      [mm10](https://www.encodeproject.org/files/ENCFF547MET/)
#'                                    only). Alternatively, a data frame or
#'                                    tibble can be provided listing the genomic
#'                                    regions to remove (having `chrom`,
#'                                    `start`, and `end`  column names). If set
#'                                    to 'NULL' (default), this step will be
#'                                    skipped (optional).
#'                                    Please note that if there are not matching
#'                                    entries in the 'chrom' columns of input
#'                                    and blacklist, an information message is
#'                                    displayed. This can happend und does not
#'                                    cause any problems with the script.
#' * `include_above_score_cutoff` -   Single numeric value that defines the
#'                                    `score` threshold above which all genomic
#'                                    regions will be retained. The `score`
#'                                    column in the peakCombiner input data
#'                                    should be non-zero for this parameter to
#'                                    be used. It is populated by
#'                                    [peakCombiner::prepare_input_regions], and
#'                                    by default takes the value of -log10(FDR)
#'                                    if possible (e.g., using a .narrowPeak
#'                                    file from MACS2 as input). Importantly,
#'                                    applying this filter retains a variable
#'                                    number of genomic regions per sample, all
#'                                    having a score greater than the
#'                                    `include_above_score_cutoff` parameter. If
#'                                    set to 'NULL' (default), this step will
#'                                    be skipped (optional).
#' * `include_top_n_scoring` -        Single numeric value that defines how many
#'                                    of the top scoring genomic regions (using
#'                                    the column `score`) are retained. All
#'                                    other genomic regions are discarded.
#'                                    Importantly, applying this filter retains
#'                                    `include_top_n_scoring` regions per
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
#' @inheritParams center_expand_regions
#'
#' @param include_by_chromosome_name
#'          * 'NULL' (default) - No chromosome name filtering will be done.
#'          * Character vector that contains chromosomes names to be retained.
#'
#' @param exclude_by_blacklist
#'          * 'NULL' (default) - No blacklist filtering will be done.
#'          * Single value for annotation of provided blacklists ('hg38' or
#'            'mm10'), which uses ENCODE blacklist to filter.
#'          * Data frame or tibble with columns `chrom`, `start`, and `end`.
#'
#' @param include_above_score_cutoff
#'          * 'NULL' (default) - No score filtering will be done.
#'          * Single numeric value that defines the `score` threshold above
#'            which all genomic regions will be retained. This results in
#'            variable number of sites per sample.
#'
#' @param include_top_n_scoring
#'          * 'NULL' (default) - No score filtering will be done.
#'          * Single numeric value representing the number of genomic regions
#'            per sample to be retained. The genomic regions are selected from
#'            highest to lowest score, and if include_top_n_scoring > number of
#'            regions, then no filtering is done.
#'
#' @param show_messages Logical value of TRUE (default) or FALSE. Defines if
#'                      info messages are displayed or not.
#'
#' @return A tibble with the columns `chrom`, `start`, `end`, `name`, `score`,
#' `strand`, `center`, `sample_name`. The definitions of these columns are
#' described in full in the [peakCombiner::prepare_input_regions] Details.
#' Use as input for functions [peakCombiner::center_expand_regions] and
#' [peakCombiner::combine_regions].
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
#' input_data <- peakCombiner::syn_data_bed
#' input_data
#'
#' data_prepared <- prepare_input_regions(
#'   data = input_data,
#'   show_messages = TRUE
#' )
#'
#' # Here use options for all four filtering methods.
#'
#' filter_regions(
#'   data = data_prepared,
#'   include_by_chromosome_name = c("chr1", "chr2", "chr4"),
#'   exclude_by_blacklist = "hg38",
#'   include_above_score_cutoff = 10,
#'   include_top_n_scoring = 100,
#'   show_messages = TRUE
#' )
#'
filter_regions <- function(data,
                           include_by_chromosome_name = NULL,
                           exclude_by_blacklist = NULL,
                           include_above_score_cutoff = NULL,
                           include_top_n_scoring = NULL,
                           show_messages = TRUE) {
  ### -----------------------------------------------------------------------###
  ### Define parameters
  ### -----------------------------------------------------------------------###
  ##
  set.seed(1234)
  ##
  ## Pass data into new variable
  data_filtered <- data

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
    data = data_filtered
  )

  ### -----------------------------------------------------------------------###
  ### Filter by chromosomes names
  ### -----------------------------------------------------------------------###

  data_filtered <-
    filter_by_chromosome_names(
      data = data_filtered,
      include_by_chromosome_name = include_by_chromosome_name
    )

  ### -----------------------------------------------------------------------###
  ### Filter by blacklist
  ### -----------------------------------------------------------------------###

  data_filtered <-
    filter_by_blacklist(
      data = data_filtered,
      exclude_by_blacklist = exclude_by_blacklist
    )

  ### -----------------------------------------------------------------------###
  ### Filter for significance values (score)
  ### -----------------------------------------------------------------------###

  data_filtered <-
    filter_by_significance(
      data = data_filtered,
      include_above_score_cutoff = include_above_score_cutoff
    )

  ### -----------------------------------------------------------------------###
  ### Select top peaks
  ### -----------------------------------------------------------------------###

  data_filtered <-
    filter_by_top_enriched(
      data = data_filtered,
      include_top_n_scoring = include_top_n_scoring
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
  ### Set message display back to default
  ### -----------------------------------------------------------------------###

  if (isFALSE(show_messages)) {
    options("rlib_message_verbosity" = "default")
  }

  return(data_filtered)
}
