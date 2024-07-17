#' Calculate expansion value from median input region size
#'
#' @description
#' Calculates the parameter `expand_by` when it was set to 'NULL' in the main
#' function. 'NULL' allows for data-driven definition of the `expand_by` value.
#' It calculates the median genomic region size of the input data and uses this
#' value like a length 1 numeric vector for expansion.
#'
#' @inheritParams center_expand_regions
#'
#' @return A vector of length 1 to define region expansion.
#'

define_expansion <- function(data = data,
                             expand_by = expand_by) {

  ### -----------------------------------------------------------------------###
  ### Pre-Check up
  ### -----------------------------------------------------------------------###

  if (!exists("expand_by")) {
    # show error message independent of parameter show_messages
    options("rlib_message_verbosity" = "default")

    cli::cli_abort(c(
      "x" = "{.arg expand_by} has to be defined."
    ))
  }

  if (!is.null(expand_by) && !is.numeric(expand_by)) {
    # show error message independent of parameter show_messages
    options("rlib_message_verbosity" = "default")

    cli::cli_abort(c(
      "x" = "{.arg expand_by} has to be {.cls numeric}."
    ))
  }

  if (any(is.na(expand_by))) {
    # show error message independent of parameter show_messages
    options("rlib_message_verbosity" = "default")

    cli::cli_abort(c(
      "x" = "{.arg expand_by} has the unallowed value {.val NA}."
    ))
  }

  if (any(expand_by < 1)) {
    # show error message independent of parameter show_messages
    options("rlib_message_verbosity" = "default")

    cli::cli_abort(c(
      "x" = "{.arg expand_by} has the unallowed value {.val 0}."
    ))
  }

  ### -----------------------------------------------------------------------###
  ### If expand_by is not provided, here we calculate the median peak size.
  ### -----------------------------------------------------------------------###

  if (any(is.null(expand_by)) == TRUE) {
    cli::cli_inform(c(
      "i" = "Input value for {.arg expand_by} is {.val NULL}. Median of all
      input genomic regions is calculated and returned for expansion."
    ))

    ## Get median peak size

    expand_by <-
      data |>
      dplyr::ungroup() |>
      dplyr::summarise(val = stats::median(.data$end - .data$start)) |>
      dplyr::mutate(val = round(.data$val / 2, 0)) |>
      dplyr::pull()

    cli::cli_inform(c(
      "v" = "{.var expand_by} was calculated from the input data and set to
      \"{.val {expand_by}}\".",
      "i" = "Genomic regions will be expanded by {expand_by}bp in
      both direction.",
      " " = " "
    ))
  } else if (is.numeric(expand_by) && length(expand_by) < 3) {
    if (length(expand_by) == 2) {
      expand_by <- round(expand_by)

      cli::cli_inform(c(
        "v" = "Parameter {.var expand_by} was defined as
        {expand_by}bp in a vector having length 2.",
        "i" = "Genomic regions will be expanded by subtracting {expand_by[1]}bp
         from, and adding {expand_by[2]}bp to, the {.field center}."
      ))
    } else if (length(expand_by) == 1) {
      expand_by <- round(expand_by)
      cli::cli_inform(c(
        "v" = "The parameter {.var expand_by} was defined as
            {expand_by}bp.",
        "i" = "Genomic regions will be expanded by {expand_by}bp
            in both directions."
      ))
    }
  } else {
    # show error message independent of parameter show_messages
    options("rlib_message_verbosity" = "default")

    cli::cli_abort(c(
      "x" = "The parameter {.var expand_by} has to be either {.val NULL} or a
          numeric vector with length 1 or 2."
    ))
  }

  ### -----------------------------------------------------------------------###
  ### Return expansion parameter
  ### -----------------------------------------------------------------------###

  return(expand_by)
}
