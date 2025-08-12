#' Calculate expansion value from median input region size
#'
#' @description
#' Calculates the parameter `expandBy` when it was set to 'NULL' in the main
#' function. 'NULL' allows for data-driven definition of the `expandBy` value.
#' It calculates the median genomic region size of the input data and uses this
#' value like a length 1 numeric vector for expansion.
#'
#' @inheritParams centerExpandRegions
#'
#' @return A vector of length 1 to define region expansion.
#'

defineExpansion <- function(data = data,
                             expandBy = expandBy) {
  ### -----------------------------------------------------------------------###
  ### Pre-Check up
  ### -----------------------------------------------------------------------###

  if (!exists("expandBy")) {
    # show error message independent of parameter showMessages
    options("rlib_message_verbosity" = "default")

    cli::cli_abort(c(
      "x" = "{.arg expandBy} has to be defined."
    ))
  }

  if (!is.null(expandBy) && !is.numeric(expandBy)) {
    # show error message independent of parameter showMessages
    options("rlib_message_verbosity" = "default")

    cli::cli_abort(c(
      "x" = "{.arg expandBy} has to be {.cls numeric}."
    ))
  }

  if (any(is.na(expandBy))) {
    # show error message independent of parameter showMessages
    options("rlib_message_verbosity" = "default")

    cli::cli_abort(c(
      "x" = "{.arg expandBy} has the unallowed value {.val NA}."
    ))
  }

  if (any(expandBy < 1)) {
    # show error message independent of parameter showMessages
    options("rlib_message_verbosity" = "default")

    cli::cli_abort(c(
      "x" = "{.arg expandBy} has the unallowed value {.val 0}."
    ))
  }

  ### -----------------------------------------------------------------------###
  ### If expandBy is not provided, here we calculate the median peak size.
  ### -----------------------------------------------------------------------###

  if (any(is.null(expandBy)) == TRUE) {
    cli::cli_inform(c(
      "i" = "Input value for {.arg expandBy} is {.val NULL}. Median of all
      input genomic regions is calculated and returned for expansion."
    ))

    ## Get median peak size

    expandBy <-
      data |>
      dplyr::ungroup() |>
      dplyr::summarise(val = stats::median(.data$end - .data$start)) |>
      dplyr::mutate(val = round(.data$val / 2, 0)) |>
      dplyr::pull()

    cli::cli_inform(c(
      "v" = "{.var expandBy} was calculated from the input data and set to
      \"{.val {expandBy}}\".",
      "i" = "Genomic regions will be expanded by {expandBy}bp in
      both direction.",
      " " = " "
    ))
  } else if (is.numeric(expandBy) && length(expandBy) < 3) {
    if (length(expandBy) == 2) {
      expandBy <- round(expandBy)

      cli::cli_inform(c(
        "v" = "Parameter {.var expandBy} was defined as
        {expandBy}bp in a vector having length 2.",
        "i" = "Genomic regions will be expanded by subtracting {expandBy[1]}bp
         from, and adding {expandBy[2]}bp to, the {.field center}."
      ))
    } else if (length(expandBy) == 1) {
      expandBy <- round(expandBy)
      cli::cli_inform(c(
        "v" = "The parameter {.var expandBy} was defined as
            {expandBy}bp.",
        "i" = "Genomic regions will be expanded by {expandBy}bp
            in both directions."
      ))
    }
  } else {
    # show error message independent of parameter showMessages
    options("rlib_message_verbosity" = "default")

    cli::cli_abort(c(
      "x" = "The parameter {.var expandBy} has to be either {.val NULL} or a
          numeric vector with length 1 or 2."
    ))
  }

  ### -----------------------------------------------------------------------###
  ### Return expansion parameter
  ### -----------------------------------------------------------------------###

  return(expandBy)
}
