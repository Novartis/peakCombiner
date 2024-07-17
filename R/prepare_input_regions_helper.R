#' Load in data based on provided in memory data frame
#'
#' @description
#' Helper function for main function [peakCombiner::prepare_input_regions].
#' Requires in memory data frame listing each sample's peak file location as
#' input and uses the provided file paths to load the input files.
#' For details see the details for [peakCombiner::prepare_input_regions].
#'
#'
#' @inheritParams prepare_input_regions
#'
#' @return A tibble with the columns `chrom`, `start`, `end`, `score`,
#' `strand`, `sample_name`. This data frame has to be further modified within
#' the function [peakCombiner::prepare_input_regions].
#'
#'
#' @noRd
load_input_regions <- function(data) {
  ### -----------------------------------------------------------------------###
  ### Define variables
  ### -----------------------------------------------------------------------###

  allowed_col_names <- c(
    "sample_name", "file_path", "file_format",
    "score_colname"
  )

  all_other_colnames <- c(
    "chrom", "start", "end", "name", "strand", "center",
    "sample_name"
  )

  score_colname <- NULL

  output_colnames <- c(
    "chrom", "start", "end", "score", "strand", "summit", "sample_name"
  )

  file_format <- data |>
    dplyr::pull(.data$file_format) |>
    unique() |>
    tolower()

  ### -----------------------------------------------------------------------###
  ### Create lookup table - HOW TO PROVIDE AS EXTERNAL DF?
  ### -----------------------------------------------------------------------###

  lookup_table_colnames_input <- tibble::tribble(
    ~file_format, ~colnames_pc,
    "narrowpeak", c(
      "chrom", "start", "end", "name", "score", "strand",
      "SignalValue", "pValue", "qValue", "peak"
    ),
    "broadpeak", c(
      "chrom", "start", "end", "name", "score", "strand",
      "SignalValue", "pValue", "qValue", "peak"
    ),
    "bed", c("chrom", "start", "end"),
    "bed6", c("chrom", "start", "end", "name", "score", "strand")
  )

  colnames_input <- lookup_table_colnames_input |>
    dplyr::filter(.data$file_format == !!file_format) |>
    dplyr::pull(.data$colnames_pc) |>
    unlist()

  ### -----------------------------------------------------------------------###
  ### Pre-Check up
  ### -----------------------------------------------------------------------###
  ## Check if data is a data_frame
  if (!is.data.frame(data)) {

    # show error message independent of parameter show_messages
    options("rlib_message_verbosity" = "default")

    cli::cli_abort(c(
      "x" = "{.arg data} must be a {.cls tbl_df/tbl/data.frame}.",
      "!" = "{.arg data} has class {.cls {class(data)}}."
    ), call. = FALSE)
  }

  ## Check if samples_sheet has correct col names
  if (!all(allowed_col_names %in% colnames(data))) {
    missing_column <- allowed_col_names[!colnames(data) %in% allowed_col_names]

    # show error message independent of parameter show_messages
    options("rlib_message_verbosity" = "default")

    cli::cli_abort(c(
      "!" = "{.arg data} must contain columns with the names
    `sample_name`, `file_path` and `file_format`.",
      "x" = "Can't find column with name '{missing_column}' in {.arg data}."
    ), call. = FALSE)

    rm(missing_column)
  }

  ### -----------------------------------------------------------------------###
  ## Test that file_format and score_colname are unique

  if (any(colnames(data) == "score_colname")) {
    n_unique_score_or_formats <-
      data |>
      dplyr::select(file_format, score_colname) |>
      unique() |>
      dplyr::count() |>
      dplyr::pull()
  } else if (!any(colnames(data) == "score_colname")) {
    n_unique_score_or_formats <-
      data |>
      dplyr::select(file_format) |>
      unique() |>
      dplyr::count() |>
      dplyr::pull()
  }

  if (n_unique_score_or_formats > 1) {

    # show error message independent of parameter show_messages
    options("rlib_message_verbosity" = "default")

    cli::cli_abort(c(
      ">" = "Argument {.arg data} can only have one unique value for
      {.field file_format} and one value for {.field score_colname}
      (if provided).",
      "x" = "Non-unique values found in {.field file_format} or
      {.field score_colname}.",
      " " = " "
    ), call. = FALSE)
  }

  remove(n_unique_score_or_formats)

  ### -----------------------------------------------------------------------###
  ## Test if provided file paths in input do exist

  if (!all(file.exists(data$file_path))) {
    # show error message independent of parameter show_messages
    options("rlib_message_verbosity" = "default")
    cli::cli_abort(c(
      ">" = "`data` contains column with name 'file_path'.",
      "x" = "At least one file does not exist."
    ), call. = FALSE)
  }

  ### -----------------------------------------------------------------------###
  ## Test if sample names are unique

  if (any(data |> dplyr::pull(.data$sample_name) |> duplicated())) {
    cli::cli_inform(c(
      "x" = "Samplesheet contains non-unique {.arg sample_name}.",
      "!" = "All {.arg sample_name} have to be unique.",
      " " = " "
    ), call. = FALSE)
  }

  ## if score_colname is not in input data, we define it
  if (!any(colnames(data) == "score_colname") ||
    any(is.null(data$score_colname)) ||
    any(is.na(data$score_colname))) {
    if (file_format %in% c("narrowpeak", "broadpeak")) {
      score_colname <- "qValue"

      cli::cli_inform(c(
        ">" = "Input data frame does not contain a column named
        {.field 'score_colname'}.",
        "i" = "Input data column {.field 'file_format'} contains the value
        {.val {file_format}}.",
        "v" = "Based on these inputs, the values from the linked regions file
        column {.field 'qValue'} will be stored in the output column
        {.field 'scores'}."
      ))
    } else if (file_format %in% c("bed")) {
      score_colname <- "score"

      cli::cli_inform(c(
        ">" = "Input data frame does not contain a column named
        {.field 'score_colname'}.",
        "i" = "Input data column {.field 'file_format'} contains the value
        {.val {file_format}}.",
        "v" = "Based on these inputs, the value {.val 0} will be stored in
        the output column {.field 'scores'}."
      ))
    }
  } else if (any(colnames(data) == "score_colname")) {
    ## if score_colname is in input data, we extract it

    score_colname_input <- data |>
      dplyr::pull(score_colname) |>
      unique()

    file_format_colnames <- lookup_table_colnames_input |>
      dplyr::filter(.data$file_format == !!file_format) |>
      dplyr::pull(.data$colnames_pc) |>
      unlist()

    if (is.numeric(score_colname_input)) {
      score_colname <- file_format_colnames[score_colname_input]
    } else if (is.character(score_colname_input)) {
      score_colname <- score_colname_input
    }
  }

  ### -----------------------------------------------------------------------###
  ### Load in peak files
  ### -----------------------------------------------------------------------###

  if (score_colname %in% all_other_colnames) {

    # show error message independent of parameter show_messages
    options("rlib_message_verbosity" = "default")

    cli::cli_abort(c(
      "x" = "Provided column name for scores {.field score_colname} is identical
      with other required output column name.",
      ">" = "Provided colnames for scores is {.val {score_colname}}.",
      ">" = "Forbidden names for {.arg score_colname} are
        {.val {required_colnames}}."
    ))
  }

  ### -----------------------------------------------------------------------###

  cli::cli_inform(c(
    ">" = "Start reading in data."
  ))
  
  ## Read in peak files
  data_readin <-
    tibble::tibble(
      sample_name = data$sample_name,
      file_path = data$file_path
    ) |>
    dplyr::mutate(
      input_file = purrr::map(.data$file_path,
        readr::read_tsv,
        col_names = colnames_input,
        .progress = FALSE,
        show_col_types = FALSE
      ) |> stats::setNames(data$sample_name)
    ) |>
    dplyr::select(-"file_path") |>
    tidyr::unnest(cols = c("input_file")) |>
    dplyr::group_by(.data$sample_name) |>
    dplyr::mutate(
      start = .data$start + 1,
      summit = .data$peak,
      score = get(score_colname)
    ) |>
    dplyr::select(paste(output_colnames, sep = ",")) |>
    dplyr::ungroup()

  ## Return the combined data
  cli::cli_inform(c(
    "v" = "Data was successfully read in."
  ))

  ### -----------------------------------------------------------------------###
  ### Return data frame
  ### -----------------------------------------------------------------------###

  return(data_readin)
}

################################################################################
################################################################################
################################################################################

#' collapse_summits
#'
#' @description
#' Helper function for main function [peakCombiner::prepare_input_regions].
#' Input data is checked for multiple entries of the same genomic
#' region. This can occur when using called peak files as multiple summits can
#' be annotated within the same genomic regions (defined by `chrom`, `start`
#' and `end`). To avoid multiple entries, this script is checking the input for
#' multiple summits within the same regions and maintains only the strongest
#' enriched (based on the values in the column `score`). This step is mandatory
#' to quantity an optimal result.
#' For details see the details for [peakCombiner::prepare_input_regions].
#'
#' @param data A tibble with the columns `chrom`, `start`, `end`, `name`,
#' `score`, `strand`, `center`, `sample_name`.
#'
#' @return A tibble with the columns `chrom`, `start`, `end`, `name`, `score`,
#' `strand`, `center`, `sample_name`. The definitions of these columns are
#' described in full in the Details below. Use as input for functions
#' [peakCombiner::center_expand_regions()], [peakCombiner::filter_regions()] and
#' [peakCombiner::combine_regions()].
#'
collapse_summits <- function(data) {
  cli::cli_inform(c(
    ">" = "Checking whether duplicated regions exist and need to be collapsed."
  ))

  ### -----------------------------------------------------------------------###
  ### Filter for top enriched regions
  ### -----------------------------------------------------------------------###
  ## Reduce to max enriched summit per region in each sample

  n_rows_before <- data$chrom |> length()

  data_prepared_filtered <- data |>
    dplyr::group_by(
      .data$sample_name,
      .data$chrom,
      .data$start,
      .data$end
    ) |>
    dplyr::filter(.data$score == .data$score |> max()) |>
    dplyr::slice_sample(n = 1) |>
    dplyr::ungroup()

  n_rows_after <- data_prepared_filtered$chrom |> length()

  ### -----------------------------------------------------------------------###
  ### Provide user feedback on collapsing summits
  ### -----------------------------------------------------------------------###

  if (n_rows_before == n_rows_after) {
    cli::cli_inform(c(
      ">" = "Checked whether duplicated regions exist and need to be
      collapsed.",
      "v" = "No regions with mutliple summits found. No regions were removed.",
      " " = " "
    ))
  } else {
    cli::cli_inform(c(
      ">" = "Checked whether duplicated regions exist and need to be
      collapsed.",
      "v" = "Duplicated regions identified and collapsed to unique
      {.field chrom}, {.field start}, and {.field end} for each sample by
      strongest {.field score} value.",
      " " = " "
    ))
  }

  ### -----------------------------------------------------------------------###
  ### Return data frame
  ### -----------------------------------------------------------------------###
  return(data_prepared_filtered)
}
