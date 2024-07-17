#' Separate genomic regions by coverage (disjoin) and filter input regions
#'
#' @description
#' Helper function for main function [peakCombiner::combine_regions].
#' Requires in memory data frame in the standard accepted format for the
#' peakCombiner package.
#' For details see the details for [peakCombiner::combine_regions].
#'
#' @details
#' Retain overlapping genomic regions that are found in at least
#' `found_in_samples` samples. In this way, you can remove rare or
#' sample-specific regions.
#'
#' @inheritParams combine_regions
#'
#' @return A tibble with the following columns: `chrom`, `start`, `end`,
#' `width`, `strand`, `revmap`, `ranking_comb_ref`, `name`, `rowname_disjoin`.
#'
cr_disjoin_filter <- function(data,
                              found_in_samples) {

  ### -----------------------------------------------------------------------###
  ### Pre-Check up
  ### -----------------------------------------------------------------------###
  ## Check if expansion exists
  if (!exists("data")) {
    # show error message independent of parameter show_messages
    options("rlib_message_verbosity" = "default")

    cli::cli_abort(c(
      "x" = "Argument {.arg data} does not exist."
    ))
  }

  if (!is.numeric(found_in_samples) ||
    is.null(found_in_samples) ||
    is.na(found_in_samples) ||
    length(found_in_samples) > 1) {
    # show error message independent of parameter show_messages
    options("rlib_message_verbosity" = "default")

    cli::cli_abort(c(
      "x" = "Argument {.arg found_in_samples} has an unaccepted value.",
      "i" = "Argument {.arg found_in_samples} has to be a positive
      {.cls numeric} value."
    ))
  }

  if (found_in_samples <= 0) {
    # show error message independent of parameter show_messages
    options("rlib_message_verbosity" = "default")

    cli::cli_abort(c(
      "x" = "Argument {.arg found_in_samples} is {.val 0} or negative.",
      "i" = "Argument {.arg found_in_samples} has to be positive."
    ))
  } else if (found_in_samples > 0 && found_in_samples < 1) {
    n_samples <-
      data |>
      dplyr::select(.data$sample_name) |>
      unique() |>
      dplyr::summarize(n = dplyr::n()) |>
      dplyr::pull(.data$n)

    found_in_samples <- ceiling(n_samples * found_in_samples)

    # TODO: Write and informational message to user about number of samples

    cli::cli_inform(c(
      ">" = "Argument {.arg found_in_samples} is between {.val 0} and {.val 1}
      and considered to be a fraction.",
      "i" = "Argument {.arg found_in_samples} was calulated based on the number
      of provided values in {.field samples_names}.",
      "v" = "Argument {.arg found_in_samples} was set to
      {.val {found_in_samples}}.",
      " " = " "
    ))

    remove(n_samples)
  }

  ### -----------------------------------------------------------------------###
  ### Combine peaks
  ### -----------------------------------------------------------------------###

  cli::cli_inform(c(
    ">" = "Start with the disjoining and filtering genomic regions."
  ))

  data_disjoin <-
    data |>
    GenomicRanges::makeGRangesFromDataFrame(
      keep.extra.columns = TRUE
    ) |>
    ## Make GRanges Object from input files -
    GenomicRanges::disjoin(with.revmap = TRUE) |>
    tibble::as_tibble() |>
    dplyr::rowwise() |>
    dplyr::mutate(count_regions = length(.data$revmap)) |>
    dplyr::ungroup() |>
    tidyr::unnest("revmap") |>
    dplyr::rename(chrom = "seqnames") |>
    dplyr::left_join(data |>
      tibble::rownames_to_column(var = "revmap") |>
      dplyr::mutate(revmap = as.integer(.data$revmap)) |>
      dplyr::select("revmap", "sample_name"), 
      by = "revmap")

  data_disjoin_meta <- data_disjoin |>
    dplyr::left_join(
      data_disjoin |>
        dplyr::group_by(.data$chrom, .data$start, .data$end) |>
        dplyr::summarise(
          length = dplyr::n_distinct(.data$sample_name),
          .groups = "drop"
        ) |>
        dplyr::select("chrom", "start", "end", "length"),
      by = c("chrom", "start", "end")
    )

  data_disjoin_filter <-
    data_disjoin_meta |>
    dplyr::filter(.data$length >= found_in_samples) |>
    tibble::rownames_to_column() |>
    dplyr::relocate(
      "ranking_comb_ref" = "rowname",
      .after = tidyr::last_col()
    ) |>
    dplyr::select(-"length", -"count_regions") |>
    tidyr::unnest("revmap", keep_empty = TRUE) |>
    dplyr::inner_join(
      data |>
        dplyr::mutate(revmap = dplyr::row_number()) |>
        dplyr::select("revmap", "name", "center", "score"),
      by = "revmap"
    ) |> ## 2: Add the metadata ID to the disjoined data frame
    tibble::rownames_to_column() |>
    dplyr::relocate(
      "rowname_disjoin" = "rowname",
      .after = tidyr::last_col()
    )

  ### -----------------------------------------------------------------------###
  ### Return prepared input data
  ### -----------------------------------------------------------------------###

  cli::cli_inform(c(
    "v" = "Disjoin and filter by {.arg found_in_samples} of genomic regions
    successfully finished.",
    " " = " "
  ))

  data_disjoin_filter <-
    data_disjoin_filter |>
    dplyr::mutate(
      start = as.numeric(.data$start),
      end = as.numeric(.data$end)
    ) |>
    dplyr::relocate("chrom", .before = "start") |>
    dplyr::ungroup()

  return(data_disjoin_filter)
}

################################################################################
################################################################################
################################################################################

#' Recombine (reduce) input regions
#'
#' @description
#' Helper function for main function [peakCombiner::combine_regions].
#' Requires in memory data frame in the standard accepted format for the
#' peakCombiner package.
#' For details see the details for [peakCombiner::combine_regions].
#'
#' @details
#' Recombine filtered genomic regions from disjoin function to create the
#' consensus regions.
#'
#'
#' @inheritParams combine_regions
#'
#' @return A tibble with the following columns: `chrom`, `start`, `end`,
#' `width`, `strand`, `name`.
#'
#'
cr_reduce <- function(data) {

  ### -----------------------------------------------------------------------###
  ### Correct parameters & load needed variables
  ### -----------------------------------------------------------------------###

  required_colnames <- c(
    "chrom", "start", "end", "width", "strand", "revmap",
    "ranking_comb_ref", "rowname_disjoin", "name"
  )

  ### -----------------------------------------------------------------------###
  ### Pre-Check up
  ### -----------------------------------------------------------------------###
  ## Check if expansion exists
  if (!exists("data")) {
    # show error message independent of parameter show_messages
    options("rlib_message_verbosity" = "default")

    cli::cli_abort(c(
      "x" = "Parameter {.arg data} does not exist."
    ))
  }

  if (!any(colnames(data) %in% required_colnames)) {
    missing_cols <- required_colnames[!required_colnames %in% colnames(data)]

    # show error message independent of parameter show_messages
    options("rlib_message_verbosity" = "default")

    cli::cli_abort(c(
      "x" = "Input dataframe {.arg data} does not contain certain columns.",
      "!" = "Missing columns might be: {missing_cols}",
      ">" = "Note to user: Please check colnames of input data.",
      ">" = "Requirde column names are: {required_colnames}"
    ))

    rm(missing_cols)
  }

  ### -----------------------------------------------------------------------###
  ### Combine peaks - Reduce
  ### -----------------------------------------------------------------------###
  ## 2: Reduce the disjoined, meta data added table
  cli::cli_inform(c(
    ">" = "Start with combining remaining genomic regions."
  ))

  data_disjoin_reduce <-
    data |>
    GenomicRanges::makeGRangesFromDataFrame(
      keep.extra.columns = TRUE,
    ) |>
    GenomicRanges::reduce(
      min.gapwidth = 1L,
      with.revmap = TRUE,
      drop.empty.ranges = TRUE
    ) |>
    tibble::as_tibble() |>
    tidyr::unnest("revmap", keep_empty = TRUE) |>
    dplyr::inner_join(
      data |>
        dplyr::mutate(revmap = dplyr::row_number()) |>
        dplyr::select("revmap", "name", "center", "score"),
      by = "revmap"
    ) |>
    dplyr::select(-"revmap") |>
    dplyr::arrange(.data$seqnames, .data$start, .data$name) |>
    dplyr::rename(chrom = .data$seqnames) |>
    unique() |>
    dplyr::ungroup()

  ### -----------------------------------------------------------------------###
  ### Return prepared input data
  ### -----------------------------------------------------------------------###
  cli::cli_inform(c(
    "v" = "Combining remaining genomic regions was successfully finished.",
    " " = " "
  ))

  return(data_disjoin_reduce)

  ### -----------------------------------------------------------------------###
}

################################################################################
################################################################################
################################################################################

#' Overlap genomic regions with original summits to remove false positive
#'
#' @description
#' Helper function for main function [peakCombiner::combine_regions].
#' Requires in memory data frame in the standard accepted format for the
#' peakCombiner package.
#' For details see the details for [peakCombiner::combine_regions].
#'
#' @details
#' Overlapping genomic regions must contain at least one 'center' from its
#' input sample regions to be considered a valid genomic region. Regions without
#' overlap might be a consequence of the expansion parameter and are likely to
#' be false positive.
#'
#' @inheritParams combine_regions
#'
#' @param input The original input file from `combine_regions` to extract center
#' information
#'
#' @return A tibble with the following columns: `chrom`, `start`, `end`,
#' `width`, `strand`, `name`.
#'
cr_overlap_with_summits <- function(data,
                                    input) {

  ### -----------------------------------------------------------------------###
  ### Correct parameters & load needed variables
  ### -----------------------------------------------------------------------###

  required_colnames <- c(
    "chrom", "start", "end", "strand", "name", "score",
    "center", "sample_name"
  )

  ### -----------------------------------------------------------------------###
  ### Pre-Check up
  ### -----------------------------------------------------------------------###
  ## Check if expansion exists
  if (!exists("input")) {
    # show error message independent of parameter show_messages
    options("rlib_message_verbosity" = "default")

    cli::cli_abort(c(
      "x" = "Argument {.arg input} does not exist."
    ))
  }

  if (!exists("data")) {
    # show error message independent of parameter show_messages
    options("rlib_message_verbosity" = "default")

    cli::cli_abort(c(
      "x" = "Argument {.arg data} does not exist."
    ))
  }

  if (!any(colnames(input) %in% required_colnames)) {
    missing_cols <- required_colnames[!required_colnames %in% colnames(input)]
    # show error message independent of parameter show_messages
    options("rlib_message_verbosity" = "default")

    cli::cli_abort(c(
      "x" = "Input dataframe {.arg input} does not contain certain columns.",
      "!" = "Missing columns might be: {missing_cols}",
      ">" = "Note to user: Please check colnames of input data.",
      ">" = "Requirde column names are: {required_colnames}"
    ))

    rm(missing_cols)
  }

  ### -----------------------------------------------------------------------###
  ### Combine peaks - Overlap with summit
  ### -----------------------------------------------------------------------###
  ## 3: Remove false positive peaks without summit
  cli::cli_inform(c(
    ">" = "Start with identification of overlaps between the original summit and
    remaining genomic regions.",
    "i" = "Remaining regions without overlap will be
    removed."
  ))

  ## Define summits from input data
  summits <-
    data |>
    dplyr::mutate(
      start = .data$center,
      end = .data$center + 1
    ) |>
    dplyr::select("chrom", "start", "end") |>
    unique()

  ## Get overlap between NEW PEAKs and original summits
  overlap_meta <-
    GenomicRanges::countOverlaps(
      GenomicRanges::makeGRangesFromDataFrame(
        data |>
          dplyr::select("chrom", "start", "end") |>
          unique(),
        keep.extra.columns = TRUE
      ),
      GenomicRanges::makeGRangesFromDataFrame(summits,
        keep.extra.columns = TRUE
      )
    ) |>
    tibble::as_tibble() |>
    tibble::rownames_to_column() |>
    dplyr::relocate("ranking" = "rowname", .after = tidyr::last_col()) |>
    dplyr::rename(count = "value")

  ## Make a df with all overlapping peaks
  data_filtered_overlap <-
    data |>
    dplyr::group_by(
      .data$chrom,
      .data$start,
      .data$end,
      .data$width,
      .data$strand
    ) |>
    tidyr::nest() |>
    tibble::rownames_to_column() |>
    dplyr::relocate("ranking" = "rowname", .after = tidyr::last_col()) |>
    dplyr::left_join(overlap_meta) |>
    dplyr::filter(!.data$count == 0) |>
    dplyr::select(-"ranking", -"count") |>
    tidyr::unnest(cols = c(data)) |>
    dplyr::ungroup()

  rm(summits, overlap_meta)

  ### -----------------------------------------------------------------------###
  ### Return prepared input data
  ### -----------------------------------------------------------------------###

  cli::cli_inform(c(
    "v" = "Retained genomic regions with input data summit overlap was
    successfully finished. ",
    " " = " "
  ))

  return(data_filtered_overlap)
}


################################################################################
################################################################################
################################################################################

#' Update the center and score information
#'
#' @description
#' Helper function for main function [peakCombiner::combine_regions].
#' Requires in memory data frame in the standard accepted format for the
#' peakCombiner package.
#' For details see the details for [peakCombiner::combine_regions].
#'
#' @details
#' As you can use the output data from this step again (e.g., to
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
#'    + 'sample_name' is a concatenated string of all input sample_names
#'
#' In addition, the output data.frame columns `sample_name`, `name` and `score`
#' will be updated.
#' 
#' @inheritParams combine_regions
#' 
#' @param input The original input file from `combine_regions` to extract center
#' information
#'
#' @return A tibble with the following columns: `chrom`, `start`, `end`, `name`,
#' `score`, `strand`, `center`, `sample_name`.
#'
cr_add_summit <- function(data,
                          input,
                          combined_center = "nearest",
                          annotate_with_input_names = FALSE,
                          combined_sample_name = NULL) {

  ### -----------------------------------------------------------------------###
  ### Correct parameters & load needed variables
  ### -----------------------------------------------------------------------###

  center_values <- c("nearest", "strongest", "middle")

  combined_center <- tolower(combined_center)

  if (!any(names(data) == "sample_name")) {
    data <- data |>
      tidyr::separate(.data$name,
        into = c("sample_name", "row_number"),
        sep = "\\|",
        remove = FALSE
      ) |>
      dplyr::select(-"row_number")
  }

  ### -----------------------------------------------------------------------###
  ### Pre-Check up
  ### -----------------------------------------------------------------------###
  ## Check if expansion exists
  if (!exists("input")) {
    # show error message independent of parameter show_messages
    options("rlib_message_verbosity" = "default")

    cli::cli_abort(c(
      "x" = "Argument {.arg input} does not exist."
    ))
  }

  if (!exists("data")) {
    # show error message independent of parameter show_messages
    options("rlib_message_verbosity" = "default")

    cli::cli_abort(c(
      "x" = "Argument {.arg data} does not exist."
    ))
  }

  ### -----------------------------------------------------------------------###

  if (!length(combined_center) == 1) {
    # show error message independent of parameter show_messages
    options("rlib_message_verbosity" = "default")

    cli::cli_abort(c(
      "x" = "Argument {.arg combined_center} is has a length of
      {length(combined_center)}.",
      ">" = "Allowed length for argument {.arg combined_center} is 1."
    ))
  }

  if (!combined_center %in% center_values) {
    # show error message independent of parameter show_messages
    options("rlib_message_verbosity" = "default")

    cli::cli_abort(c(
      "x" = "Argument {.arg combined_center} is an unallowed value:
      {.val {combined_center}}.",
      ">" = "Allowed values for argument {.arg combined_center} are:
      {.val {center_values}}."
    ))
  }

  ### -----------------------------------------------------------------------###
  if (!is.logical(annotate_with_input_names)) {
    # show error message independent of parameter show_messages
    options("rlib_message_verbosity" = "default")

    cli::cli_abort(c(
      "x" = "Argument {.arg annotate_with_input_names} has to be
    {.cls logical}.",
      "i" = "Accepted values are either {.val TRUE} or {.val FALSE}."
    ))
  } else if (!length(annotate_with_input_names) == 1) {
    # show error message independent of parameter show_messages
    options("rlib_message_verbosity" = "default")

    cli::cli_abort(c(
      "x" = "Argument {.arg annotate_with_input_names} has to be a single
    {.cls logical} value.",
      "i" = "Accepted values are either {.val TRUE} or {.val FALSE}."
    ))
  } else if (is.na(annotate_with_input_names) ||
    is.null(annotate_with_input_names)) {
    # show error message independent of parameter show_messages
    options("rlib_message_verbosity" = "default")

    cli::cli_abort(c(
      "x" = "Argument {.arg annotate_with_input_names} is an unallowed value:
      {.val {annotate_with_input_names}}.",
      "i" = "Argument {.arg annotate_with_input_names} has to be logical.
      Accepted values are either {.val TRUE} or {.val FALSE}."
    ))
  }

  ### -----------------------------------------------------------------------###

  if (is.data.frame(combined_sample_name) ||
    tibble::is_tibble(combined_sample_name)) {
    # show error message independent of parameter show_messages
    options("rlib_message_verbosity" = "default")

    cli::cli_abort(c(
      "x" = "Argument {.arg combined_sample_name} is an unallowed class:
      {.cls {class(combined_sample_name)}}.",
      ">" = "Argument {.arg combined_sample_name} has be either
      {.cls NULL} or {.cls character} value of length 1."
    ))
  } else if (!is.null(combined_sample_name) &&
    !is.character(combined_sample_name)) {
    # show error message independent of parameter show_messages
    options("rlib_message_verbosity" = "default")

    cli::cli_abort(c(
      "x" = "Argument {.arg combined_sample_name} is an unallowed value:
      {.val {combined_sample_name}}.",
      ">" = "Argument {.arg combined_sample_name} has be either {.cls NULL} or
      single {.cls character} value."
    ))
  } else if (is.character(combined_sample_name) &&
    !length(combined_sample_name) == 1) {
    # show error message independent of parameter show_messages
    options("rlib_message_verbosity" = "default")

    cli::cli_abort(c(
      "x" = "Argument {.arg combined_sample_name} is an unallowed value with
      length other then 1: {.val {combined_sample_name}}.",
      ">" = "Argument {.arg combined_sample_name} has be either {.cls NULL} or
      single {.cls character} value."
    ))
  }

  ### -----------------------------------------------------------------------###
  ## If mean of all centers is 0. combined_center has to be middle
  if (mean(data$center) == 0) {
    cli::cli_inform(c(
      "!" = "Column {.field center} in input data contains only values of
      {.val 0}.",
      ">" = "Argument {.arg combined_center} is set to {.val middle} as no
      useful center information can be identified input data.",
      " " = " "
    ))

    combined_center <- "middle"
  }

  ### -----------------------------------------------------------------------###
  ##
  cli::cli_inform(c(
    ">" = "Information from input {.field center} will be added to output
    data frame."
  ))

  data_center_add <- data |>
    dplyr::group_by(.data$chrom, .data$start, .data$end) |>
    dplyr::mutate(peak = paste("peak", seq_len(dplyr::n()), sep = "_")) |>
    dplyr::ungroup()

  if (combined_center == "nearest") {
    cli::cli_inform(c(
      "i" = "Argument {.arg combined_center} was defined as {.val nearest}.",
      "i" = "The mean of all input centers is calculated and the
    nearest input {.field center} is used",
      ">" = "Center information in {.field center} and {.field score} are added
      to the output data frame."
    ))

    data_center_add <- data_center_add |>
      dplyr::group_by(.data$chrom, .data$start, .data$end) |>
      dplyr::mutate(
        mean_center = .data$center |> mean() |> round(0),
        distance = (.data$mean_center - .data$center) |> abs()
      ) |>
      dplyr::slice_min(n = 1, order_by = .data$distance, with_ties = TRUE) |>
      dplyr::slice_max(n = 1, order_by = .data$score, with_ties = FALSE) |>
      dplyr::sample_n(size = 1) |>
      dplyr::ungroup() |>
      dplyr::select(-"width", -"peak", -"mean_center", -"distance")
  } else if (combined_center == "strongest") {
    cli::cli_inform(c(
      "i" = "Argument {.arg combined_center} was defined as {.val strongest}.",
      "i" = "Based on column {.field score} the strongest input {.field center}
      is idenfied.",
      ">" = "Center information in {.field center} and {.field score} are added
      to the output data frame."
    ))

    data_center_add <- data_center_add |>
      dplyr::group_by(.data$chrom, .data$start, .data$end) |>
      dplyr::slice_max(n = 1, order_by = .data$score, with_ties = FALSE) |>
      dplyr::sample_n(size = 1) |>
      dplyr::ungroup() |>
      dplyr::select(-"width", -"peak")
  } else if (combined_center == "middle") {
    cli::cli_inform(c(
      "i" = "Argument {.arg combined_center} was defined as {.val middle}.",
      "i" = "The middle between {.field start} and {.field end} and the mean
      {.field score} is calculated.",
      ">" = "Newly calculated center information are added to {.field center}
      and {.field score} in the output data frame."
    ))

    data_center_add <- data_center_add |>
      dplyr::group_by(.data$chrom, .data$start, .data$end) |>
      dplyr::mutate(
        center = .data$start + ((.data$end - .data$start) / 2),
        score = mean(.data$score),
        name = NA,
        sample_name = NA
      ) |>
      dplyr::ungroup() |>
      dplyr::select(-"width", -"peak") |>
      unique()
  }

  cli::cli_inform(c(
    "v" = "Output data frame columns {.field center} and {.field score} were
  successfully populated.",
    " " = " "
  ))

  ### -----------------------------------------------------------------------###

  if (is.null(combined_sample_name)) {
    cli::cli_inform(c(
      ">" = "No value for column {.field sample_name} was provided.",
      "i" = "Column {.field sample_name} is filled with all input
      sample_names.",
      "i" = "Column {.field name} is created as unique identifier for each row
      containing {.field sample_name} and the row number."
    ))

    data_center_add_out <- data_center_add |>
      dplyr::mutate(
        sample_name = input |>
          dplyr::pull(.data$sample_name) |>
          unique() |>
          paste(collapse = "|"),
        name = paste(.data$sample_name, dplyr::row_number(), sep = "|")
      ) |>
      dplyr::relocate("score", .after = "name") |>
      dplyr::relocate("sample_name", .after = "center")
  } else if (is.character(combined_sample_name) &&
    length(combined_sample_name) == 1) {
    cli::cli_inform(c(
      ">" = "The value {.val {combined_sample_name}} for column
    {.field sample_name} was provided.",
      "i" = "Column {.field sample_name} is filled with provided value
    {combined_sample_name}.",
      "i" = "Column {.field name} is created as unique identifier for each row
    containing {.field sample_name} and the row number."
    ))

    data_center_add_out <- data_center_add |>
      dplyr::mutate(
        sample_name = combined_sample_name,
        name = paste(.data$sample_name, dplyr::row_number(), sep = "|")
      ) |>
      dplyr::relocate("name", .before = "score") |>
      dplyr::relocate("center", .before = "sample_name")
  } else {
    # show error message independent of parameter show_messages
    options("rlib_message_verbosity" = "default")

    cli::cli_abort(c(
      ">" = "The value '{combined_sample_name}' for column
    {.field sample_name} was provided but is a non-accepted values."
    ))
  }

  cli::cli_inform(c(
    "v" = "The columns {.field sample_name} and {.field name} were successfully
  populated.",
    " " = " "
  ))

  ### -----------------------------------------------------------------------###

  if (isTRUE(annotate_with_input_names)) {
    cli::cli_inform(c(
      "i" = "Argument {.arg annotate_with_input_names} was set to {.val TRUE}.",
      ">" = "Column {.field input_names} is added to output data frame."
    ))

    data_input_names <- data |>
      dplyr::select("chrom", "start", "end", "width", "strand", "name") |>
      dplyr::group_by(
        .data$chrom,
        .data$start,
        .data$end,
        .data$width,
        .data$strand
      ) |>
      tidyr::nest()

    data_center_add_out <- data_center_add_out |>
      dplyr::left_join(data_input_names,
        by = c("chrom", "start", "end", "strand")
      ) |>
      dplyr::mutate(
        input_names = purrr::map_vec(data, ~ .x |>
          dplyr::pull(name) |>
          paste0(collapse = ";"))
      ) |>
      dplyr::select(
        "chrom", "start", "end", "strand", "name",
        "score", "center", "sample_name", "input_names"
      ) |>
      dplyr::mutate(
        chrom = as.character(.data$chrom),
        start = as.numeric(.data$start),
        end = as.numeric(.data$end)
      ) |>
      dplyr::relocate("input_names", .after = "sample_name")

    cli::cli_inform(c(
      "v" = "Additional column {.field input_names} was successfully
      populated.",
      " " = " "
    ))
  } else if (isFALSE(annotate_with_input_names)) {
    cli::cli_inform(c(
      "i" = "Argument {.arg annotate_with_input_names} was set to
      {.val FALSE}.",
      "v" = "Column {.field input_names} is not added to output data frame.",
      " " = " "
    ))

    # Clean up and check for class
    data_center_add_out <- data_center_add_out |>
      dplyr::select(
        "chrom", "start", "end", "strand", "name",
        "score", "center", "sample_name"
      ) |>
      dplyr::mutate(
        chrom = as.character(.data$chrom),
        start = as.numeric(.data$start),
        end = as.numeric(.data$end)
      )
  }

  ### -----------------------------------------------------------------------###

  data_center_add_out |>
    dplyr::mutate(strand = ifelse(.data$strand == "*", ".", .data$strand))

  cli::cli_inform(c(
    "v" = "All required columns in output data frame were successfully
  populated.",
    " " = " "
  ))

  return(data_center_add_out)
}
