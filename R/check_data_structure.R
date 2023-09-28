#' Control structure of peakCombiner data structure
#'
#' @description
#' This is a general helper function for the package [peakCombiner]. Aim of this
#' function is to check a data frame for the correct column names and classes of
#' each column to ensure to be an accepte inpuut for functions:
#' [peakCombiner::center_expand_regions()], [peakCombiner::filter_regions()] and
#' [peakCombiner::combine_regions()].
#'
#'
#' @param data A tibble with the columns `chrom`, `start`, `end`, `name`,
#' `score`, `strand`, `center`, `sample_name`. Additional columns are tolerated.
#'
#' @return A tibble with the columns `chrom`, `start`, `end`, `name`, `score`,
#' `strand`, `center`, `sample_name`. The definitions of these columns are
#' described in full in the Details below. Use as input for functions
#' [peakCombiner::center_expand_regions()], [peakCombiner::filter_regions()] and
#' [peakCombiner::combine_regions()].
#'
check_data_structure <- function(data) {

  ### -----------------------------------------------------------------------###
  ### Define variables
  ### -----------------------------------------------------------------------###

  accepted_strand_values <- c("+", "-", ".")

  ### -----------------------------------------------------------------------###
  ### Check if required input parameters were provided
  ### -----------------------------------------------------------------------###

  if (!exists("data")) {
    # show error message independent of parameter show_messages
    options("rlib_message_verbosity" = "default")

    cli::cli_abort(c(
      "x" = "{.arg data} does not exist."
    ))
  }

  ### -----------------------------------------------------------------------###
  ### Check if any columns is NA
  ### -----------------------------------------------------------------------###

  cli::cli_inform(c(
    ">" = "Checking {.cls class} and {.val values} of all columns."
  ))

  cols_w_na <- which(
    data |>
      is.na() |>
      colSums()
    > 0
  )

  if (length(cols_w_na) > 0) {
    # show error message independent of parameter show_messages
    options("rlib_message_verbosity" = "default")

    cli::cli_abort(c(
      "x" = "Data contains {.val NA} values in {length(cols_w_na)} columns.",
      "!" = "The following column{?s} contain{?s/} {.val NA}{?s}:
      {.field {names(cols_w_na)}}",
      ">" = "{.emph Note}: Please check data and remove {.val NA} values."
    ))
  }

  ### -----------------------------------------------------------------------###
  ### Check structure of accepted data frame
  ### -----------------------------------------------------------------------###

  ## Check chrom

  if (!is.character(data$chrom)) {
    # show error message independent of parameter show_messages
    options("rlib_message_verbosity" = "default")

    cli::cli_alert(c(
      ">" = "Column {.field 'chrom'} is not class {.cls character}. It will be
      converted to class {.cls character}."
    ))

    data <-
      data |>
      dplyr::mutate(chrom = as.character(.data$chrom))
  }

  ### -----------------------------------------------------------------------###
  ## Check start

  if (!is.numeric(data$start)) {
    # show error message independent of parameter show_messages
    options("rlib_message_verbosity" = "default")

    cli::cli_alert(c(
      ">" = "Column {.field 'start'} is not class {.cls numeric}. It will be
      converted to class {.cls numeric}."
    ))

    data <-
      data |>
      dplyr::mutate(start = as.numeric(.data$start))
  }

  ### -----------------------------------------------------------------------###
  ## Check end

  if (!is.numeric(data$end)) {
    # show error message independent of parameter show_messages
    options("rlib_message_verbosity" = "default")

    cli::cli_alert(c(
      "x" = "Column {.field 'end'} is not {.cls numeric}. It will be converted
      to {.cls numeric}."
    ))

    data <-
      data |>
      dplyr::mutate(end = as.numeric(.data$end))
  }

  ## Check if any end is before start coordinate

  n_negative_widths <-
    data |>
    dplyr::mutate(width = .data$end - .data$start) |>
    dplyr::filter(.data$width < 0) |>
    dplyr::count() |>
    dplyr::pull()

  if (n_negative_widths > 0) {
    data <-
      data |>
      dplyr::filter((.data$end - .data$start) >= 0)

    # show error message independent of parameter show_messages
    options("rlib_message_verbosity" = "default")

    cli::cli_alert(c(
      "x" = "Column {.field 'end'} is smaller than column 'start' in
      {n_negative_widths} rows. These rows have been deleted."
    ))
  }

  remove(n_negative_widths)

  ### -----------------------------------------------------------------------###
  ## Check name

  if (!is.character(data$name)) {
    # show error message independent of parameter show_messages
    options("rlib_message_verbosity" = "default")

    cli::cli_alert(c(
      ">" = "Column {.field 'name'} is not class {.cls name}. It will be
      converted to class {.cls name}."
    ))

    data <-
      data |>
      dplyr::mutate(name = as.character(.data$name))
  }

  ### -----------------------------------------------------------------------###
  ## Check for score class
  if (!is.numeric(data$score)) {
    # show error message independent of parameter show_messages
    options("rlib_message_verbosity" = "default")

    cli::cli_alert(c(
      "x" = "Column {.field 'score'} is not class {.cls numeric}. It will be
      converted to class {.cls numeric}."
    ))
    data <-
      data |>
      dplyr::mutate(score = as.numeric(.data$score))
  }

  ### -----------------------------------------------------------------------###
  ## Check for strand class
  if (!is.character(data$strand)) {
    # show error message independent of parameter show_messages
    options("rlib_message_verbosity" = "default")

    cli::cli_alert(c(
      "x" = "Column {.field 'strand'} is not class {.cls character}. It will be
      converted to class {.cls character}."
    ))

    data <-
      data |>
      dplyr::mutate(strand = as.character(.data$strand))
  }

  ## Count how often a not accept value is found in column strand
  n_unaccepted_strand_values <-
    data |>
    dplyr::filter(!.data$strand %in% accepted_strand_values) |>
    dplyr::count() |>
    dplyr::pull()

  if (n_unaccepted_strand_values > 0) {
    # show error message independent of parameter show_messages
    options("rlib_message_verbosity" = "default")

    cli::cli_abort(c(
      "x" = "Column {.field strand} contains {n_unaccepted_strand_values}
      row{?s} with unexpected values.",
      ">" = "Acceptable values for strand are {.val {accepted_strand_values}}."
    ))

    data <-
      data |>
      dplyr::mutate(strand = ifelse(is.na(.data$strand), ".", .data$strand))
  }

  remove(n_unaccepted_strand_values)

  ### -----------------------------------------------------------------------###
  ## Check for center class
  if (!is.numeric(data$center)) {
    cli::cli_inform(c(
      "!" = "Column {.field 'center'} is not of class {.cls numeric}.",
      ">" = "It will be converted to class {.cls numeric}."
    ))

    data <-
      data |>
      dplyr::mutate(center = as.numeric(.data$center))
  }

  ### -----------------------------------------------------------------------###
  ## Check for sample_name class
  if (!is.character(data$sample_name)) {
    cli::cli_inform(c(
      "!" = "Column {.field 'sample_name'} is not of class {.cls character}.",
      ">" = "It will be converted to class {.cls character}."
    ))

    data <-
      data |>
      dplyr::mutate(sample_name = as.character(.data$sample_name))
  }
  ### -----------------------------------------------------------------------###
  ## Return data

  cli::cli_inform(c(
    "v" = "Structure of data was successfully checked to be an accepted input.",
    " " = " "
  ))

  return(data)
}
