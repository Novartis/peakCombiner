#' Define which chromosomes are included
#'
#' @description
#' Retains only chromosomes that are in the provided vector. By not including
#' mitochondrial, sex, or non-classical chromosomes, genomic regions found on
#' these chromosomes can be removed. If set to 'NULL' (default), this step will
#' be skipped.
#'
#' @inheritParams filter_regions
#'
#' @return Data frame filtered by chromosome names based on the provided
#' parameters
#'
#' @noRd
#'
filter_by_chromosome_names <- function(data,
                                       include_by_chromosome_name = NULL) {
  ### -----------------------------------------------------------------------###
  ### Pre-Check up
  ### -----------------------------------------------------------------------###
  ##
  set.seed(1234)
  ##
  ## check if input vector is numeric and if so change to character
  if (!is.character(include_by_chromosome_name) &&
    !is.null(include_by_chromosome_name)) {
    cli::cli_inform(c(
      "!" = "{.arg include_by_chromosome_name} has the wrong class.",
      ">" = "{.arg include_by_chromosome_name} is converted to
      {.cls character}."
    ))

    include_by_chromosome_name <- as.character(include_by_chromosome_name)
  }

  if (!is.null(include_by_chromosome_name) &&
    !is.vector(include_by_chromosome_name)) {
    # show error message independent of parameter show_messages
    options("rlib_message_verbosity" = "default")

    cli::cli_abort(c(
      "x" = "{.arg include_by_chromosome_name} has to be a
        {.cls character} vector or {.val NULL}.",
      "!" = "Provided dataset is a
        {.cls {class(include_by_chromosome_name)}}."
    ))
  } else if (length(include_by_chromosome_name) > 1 &&
    any(is.null(include_by_chromosome_name))) {
    # show error message independent of parameter show_messages
    options("rlib_message_verbosity" = "default")

    cli::cli_abort(c(
      "x" = "{.arg include_by_chromosome_name} has a length of
      {length(include_by_chromosome_name)} and contains {.val NULL}.",
      "i" = " Allowed values for {.arg include_by_chromosome_name} are
      either a single {.val NULL} (length = 1) or vector with chromosome names
      (length >= 1)."
    ))
  } else if (any(is.na(include_by_chromosome_name))) {
    # show error message independent of parameter show_messages
    options("rlib_message_verbosity" = "default")

    cli::cli_abort(c(
      "x" = "{.arg include_by_chromosome_name} has a length of
      {length(include_by_chromosome_name)} and contains {.val NA}.",
      "i" = " Allowed values for {.arg include_by_chromosome_name} are
      either a single {.val NA} (length = 1) or vector with chromosome names
      (length >= 1)."
    ))
  } else if (is.null(include_by_chromosome_name)) {
    cli::cli_inform(c(
      "i" = "The argument {.arg include_by_chromosome_name} is {.val NULL}.",
      "v" = "No filtering for chromosome names in {.field chrom} is done.",
      " " = " "
    ))

    return(data)
  } else if (is.vector(include_by_chromosome_name)) {
    # looks good so let's move on and check the goodness of the input vector

    cli::cli_inform(c(
      ">" = "The argument {.arg include_by_chromosome_name} is a class
      {.cls character} of length {length(include_by_chromosome_name)} and will
      be used to retain matchhing chromsome names in {.field chrom}."
    ))

    # get the chromosome names that are used
    data_chr <-
      data |>
      dplyr::pull(.data$chrom) |>
      unique()

    # get the chromosome names from the user that are NOT found
    # as a chromosome name in the data
    not_found_chr_names <-
      include_by_chromosome_name[!unique(include_by_chromosome_name) %in%
        data_chr] |>
      unique()

    if (length(not_found_chr_names) > 0) {
      cli::cli_inform(c(
        "!" = "Input for {.arg include_by_chromosome_name} contains values not
            found in the input data.",
        "i" = "The following chromosome name{?s} you entered {?is/are}
               not used: {.val {not_found_chr_names}}.",
        ">" = "{.emph Is {.val include_by_chromosome_name} correctly defined?}"
      ))
    }
    rm(not_found_chr_names)

    # Get the chromosomes that were will be filtered out to report
    not_retained_chr_names <-
      data_chr[!data_chr %in% include_by_chromosome_name] |> unique()

    retained_chr_names <-
      data_chr[data_chr %in% include_by_chromosome_name] |> unique()

    ## Here is where the filtering happens!
    data <- data |>
      dplyr::filter(.data$chrom %in% retained_chr_names) |>
      dplyr::ungroup()

    cli::cli_inform(c(
      "v" = "Entries in {.field chrom} with the value{?s}
      {.val {retained_chr_names}} are retained."
    ))

    if (length(not_retained_chr_names) > 0) {
      cli::cli_inform(c(
        "i" = "The following {length(not_retained_chr_names)} entr{?y/ies} in
        column {.field chrom} from the input data {?was/were} not retained:
        {.val {not_retained_chr_names}}."
      ))
    }

    cli::cli_inform(c(
      "v" = "Input data was filtered to retain regions on defined chromosome.",
      " " = " "
    ))

    rm(not_retained_chr_names)
  }

  return(data)
}

################################################################################
################################################################################
################################################################################

#' Define if and which blacklisted regions are excluded
#'
#' @description
#' Remove ENCODE-annotated blacklisted regions for either human
#' ([hg38][https://www.encodeproject.org/files/ENCFF356LFX/]) or mouse
#' ([mm10][https://www.encodeproject.org/files/ENCFF547MET/]). Alternatively, a
#' data frame or tibble can be provided listing the genomic regions to remove
#' (having `chrom`, `start`, and `end`  column names). If set to 'NULL'
#' (default), this step will be skipped.
#'
#' @inheritParams filter_regions
#'
#' @return Data frame filtered by blacklist based on the provided parameters.
#'
#' @noRd
#'
filter_by_blacklist <- function(data,
                                exclude_by_blacklist = NULL) {
  ### -----------------------------------------------------------------------###
  ### Define parameters
  ### -----------------------------------------------------------------------###
  ##
  set.seed(1234)
  ##
  
  allowed_blacklist_annotations <- c("hg38", "mm10")
  required_colnames_blacklist <- c("chrom", "start", "end")

  ### -----------------------------------------------------------------------###
  ### Pre-Check up
  ### -----------------------------------------------------------------------###

  if (!exists("exclude_by_blacklist")) {
    # show error message independent of parameter show_messages
    options("rlib_message_verbosity" = "default")

    cli::cli_abort(c(
      "x" = "Parameter {.arg exclude_by_blacklist} doesn't exists.",
      "i" = "Allowed values are {.val {c('NULL', 'hg38', 'mm10')}} or
      a data frame with genomic coordinates in columns named {.field {c('chrom',
      'start', 'end')}}."
    ))
  }

  ### -----------------------------------------------------------------------###
  if (is.null(exclude_by_blacklist)) {
    cli::cli_inform(c(
      "i" = "The argument {.arg exclude_by_blacklist} is {.val NULL}.",
      "v" = "No filtering by blacklisted regions is done.",
      " " = " "
    ))

    return(data)
  } else if (is.data.frame(exclude_by_blacklist)) {
    ## Check for correct colnames
    colnames(exclude_by_blacklist) <- tolower(colnames(exclude_by_blacklist))

    ## Check if data frame has chrom, start and end
    if (!all(required_colnames_blacklist %in% names(exclude_by_blacklist))) {
      # show error message independent of parameter show_messages
      options("rlib_message_verbosity" = "default")

      cli::cli_abort(c(
        "x" = "{.arg exclude_by_blacklist} is a data frame and misses
          some required column names.",
        "i" = "Required column names are:
        {.field {required_colnames_blacklist}}."
      ))
    }
    ## Check if chrom is character, start and end are numeric
    if (!is.character(exclude_by_blacklist$chrom) &&
      !is.numeric(exclude_by_blacklist$start) &&
      !is.numeric(exclude_by_blacklist$end)) {
      # show error message independent of parameter show_messages
      options("rlib_message_verbosity" = "default")

      cli::cli_abort(c(
        "x" = "{.arg exclude_by_blacklist} is a data frame and some columns
          are the wrong data type",
        "!" = "Required types for columns are: {.field chrom} as
        {.cls character}, {.field start} as {.cls numeric}, {.field end} as
        {.cls numeric} ."
      ))
    }

    cli::cli_inform(c(
      ">" = "User provied dataframe will be used for blacklist filtering."
    ))

    blacklist_data <- exclude_by_blacklist
  } else if (is.character(exclude_by_blacklist) &&
    length(exclude_by_blacklist) == 1) {
    exclude_by_blacklist <- tolower(exclude_by_blacklist)

    # Check if the value is allowed and abort if not
    if (!exclude_by_blacklist %in% allowed_blacklist_annotations) {
      # show error message independent of parameter show_messages
      options("rlib_message_verbosity" = "default")

      cli::cli_abort(c(
        "x" = "{.arg exclude_by_blacklist} value for vector of length 1
        is not allowed.",
        ">" = "{.arg exclude_by_blacklist} has the value
        {.val {exclude_by_blacklist}}'.",
        "i" = "Allowed values are {.val {c('hg38', 'mm10')}}."
      ))
    }

    cli::cli_inform(c(
      ">" = "Blacklist for annotation {.val {exclude_by_blacklist}} will be
      used for filtering."
    ))

    # Load the blacklist corresponding to the character parameter hg38 or mm10
    blacklist_data <- if (exclude_by_blacklist == "hg38") {
      blacklist_hg38  <- NULL
      data("blacklist_hg38", package = "peakCombiner", envir = environment())
      blacklist_hg38
    } else if (exclude_by_blacklist == "mm10") {
      blacklist_mm10  <- NULL
      data("blacklist_mm10", package = "peakCombiner", envir = environment())
      blacklist_mm10
    } else {
      stop("Invalid genome parameter. Please use 'hg38' or 'mm10'.")
    }
    
  } else {
    # show error message independent of parameter show_messages
    options("rlib_message_verbosity" = "default")

    cli::cli_abort(c(
      "x" = "Parameter {.arg exclude_by_blacklist} is not an accepted format",
      "i" = "Allowed values are {.val {c('NULL', 'hg38', 'mm10')}} or
      a data frame with genomic coordinates."
    ))
  }

  ### -----------------------------------------------------------------------###
  ### Filter by blacklist
  ### -----------------------------------------------------------------------###

  ## Note: blacklist_data is either provided by the user or loaded from the file

  ## Get chrom names from blacklist
  blacklist_chr <-
    blacklist_data |>
    dplyr::pull(.data$chrom) |>
    unique()

  ## Get chrom names from input data
  data_chr <-
    data |>
    dplyr::pull(.data$chrom) |>
    unique()

  not_found_blacklist <- dplyr::setdiff(data_chr, blacklist_chr)
  not_found_input <- dplyr::setdiff(blacklist_chr, data_chr)

  if (length(not_found_blacklist) > 0) {
    cli::cli_inform(c(
      "!" = "Provided blacklist contains chromosome names (in {.field chrom})
    not found in input data.",
      "i" = "The following blacklist chromosome{?s} {?has/have} no match:
    {.val {not_found_blacklist}}.",
      ">" = "{.emph Note to user: Please doublecheck this observation.}"
    ))
  }

  if (length(not_found_input) > 0) {
    cli::cli_inform(c(
      "!" = "Provided input data contains chromosome names not found
    in blacklist.",
      "i" = "The following input chromosome{?s} {?has/have} no match:
    {not_found_input}.",
      ">" = "{.emph Note to user: Please doublecheck this observation.}"
    ))
  }
  rm(not_found_blacklist, not_found_input)


  ## Do filtering: match of CHR between input and blacklist
  data <-
    data |>
    GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE) |>
    IRanges::subsetByOverlaps(
      blacklist_data |>
        GenomicRanges::makeGRangesFromDataFrame(
          keep.extra.columns = TRUE, 
        ),
      invert = TRUE
    ) |> 
    suppressWarnings() |> #Recently added to solve warning
    tibble::as_tibble() |>
    dplyr::rename(chrom = "seqnames") |>
    dplyr::select(-"width") |>
    dplyr::mutate(
      chrom = as.character(.data$chrom),
      strand = as.character(.data$strand)
    ) |>
    dplyr::ungroup()

  cli::cli_inform(c(
    "v" = "Input data was filtered by blacklist.",
    " " = " "
  ))


  ### -----------------------------------------------------------------------###
  ### Return data frame
  ### -----------------------------------------------------------------------###

  return(data)
}

################################################################################
################################################################################
################################################################################

#' Define a score/significance cutoff above which regions are included
#'
#' @description
#' Single numeric value that defines the `score` threshold above which all
#' genomic regions will be retained. The `score` column in the peakCombiner
#' input data should be non-zero for this parameter to be used. It is populated
#' by [peakCombiner::prepare_input_regions], and by default takes the value of
#' -log10(FDR) if possible (e.g., using a .narrowPeak file from MACS2 as input).
#' Importantly, applying this filter retains a variable number of genomic
#' regions per sample, all having a score greater than the
#' `include_above_score_cutoff` parameter. If set to 'NULL' (default), this
#' step will be skipped.
#'
#' @inheritParams filter_regions
#'
#' @noRd
#'
filter_by_significance <- function(data,
                                   include_above_score_cutoff = NULL) {
  ##
  set.seed(1234)
  ##
  if (is.null(include_above_score_cutoff)) {
    cli::cli_inform(c(
      "i" = "The argument {.arg include_above_score_cutoff} is {.val NULL}.",
      "v" = "No filtering by {.field score} threshold is done.",
      " " = " "
    ))

    return(data)
  } else if (is.numeric(include_above_score_cutoff) &&
    length(include_above_score_cutoff) == 1) {
    cli::cli_inform(c(
      ">" = "Significance in {.field score} is filtered and all regions above
      {.val {include_above_score_cutoff}} will be retained."
    ))

    ## Format looks good, filter by the include_above_score_cutoff value

    input_rows <- nrow(data)

    data <-
      data |>
      dplyr::arrange(.data$sample_name, .data$score) |>
      dplyr::filter(.data$score >= include_above_score_cutoff) |>
      dplyr::ungroup()

    cli::cli_inform(c(
      "i" = "A total of {nrow(data)} of {input_rows} input regions are
      retained with value in {.field score} a above
      {include_above_score_cutoff}. ",
      "v" = "Input data was filtered to retain regions with a {.field score}
      above the defined threshold.",
      " " = " "
    ))
    rm(input_rows)
  } else {
    # show error message independent of parameter show_messages
    options("rlib_message_verbosity" = "default")

    cli::cli_abort(c(
      "x" = "Argument 'include_above_score_cutoff' has to be a
        {.cls numeric}.",
      "!" = "Provided argument is a
        {.cls {class(include_above_score_cutoff)}}"
    ))
  }

  ### -----------------------------------------------------------------------###
  ### Return data
  ### -----------------------------------------------------------------------###

  return(data)
}

################################################################################
################################################################################
################################################################################

#' Define if a certain number of regions per sample are included
#'
#' @description
#' Single numeric value that defines how many of the top scoring genomic regions
#' (using the column `score`) are retained. All other genomic regions are
#' discarded. Importantly, applying this filter retains `include_top_n_scoring`
#' regions per sample, which means that the minimum enrichment levels may vary
#' between samples. Note that if multiple genomic regions have the same `score`
#' cutoff value, then all of those genomic regions are included. In this case,
#' the number of resulting regions retained may be a bit higher than the input
#' parameter. If set to 'NULL' (default), this step will be skipped.
#'
#' @inheritParams filter_regions
#'
#' @return Data frame filtered by top enriched regions based on the
#' provided parameters.
#'
#' @noRd
#'
filter_by_top_enriched <- function(data,
                                   include_top_n_scoring = include_top_n_scoring) {
  ##
  set.seed(1234)
  ##
  if (is.null(include_top_n_scoring)) {
    cli::cli_inform(c(
      "i" = "The argument {.arg include_top_n_scoring} is {.val NULL}.",
      "v" = "No top enriched regions were selected. All input regions are
      retained.",
      " " = " "
    ))

    return(data)
  } else if (is.numeric(include_top_n_scoring) &&
    include_top_n_scoring > 0) {
    ### ---------------------------------------------------------------------###

    cli::cli_inform(c(
      "i" = "The argument {.arg include_top_n_scoring} extracted the the top
      {.num {include_top_n_scoring}} regions by {.field score} per sample (based
      on the values in {.field sample_name}).",
      ">" = "The top enriched {.num {include_top_n_scoring}} regions per sample
      will be retained."
    ))

    too_few_regions_left <-
      data |>
      dplyr::group_by(.data$sample_name) |>
      dplyr::summarise(counts = dplyr::n(), .groups = "drop") |>
      dplyr::filter(.data$counts < include_top_n_scoring) |>
      dplyr::pull(.data$sample_name)

    if (length(too_few_regions_left) > 0) {
      cli::cli_inform(c(
        "i" = "The argument {.arg include_top_n_scoring} was defined as
        {include_top_n_scoring}.",
        ">" = "The following {.val sample_names} contain less regions then
        defined by {.arg include_top_n_scoring}: {too_few_regions_left}",
        "!" = "No genomic regions will be removed for such samples."
      ))
    }

    ### ---------------------------------------------------------------------###

    data <-
      data |>
      dplyr::group_by(.data$sample_name) |>
      dplyr::top_n(n = !!include_top_n_scoring, wt = .data$score) |>
      dplyr::ungroup()

    cli::cli_inform(c(
      "v" = "Input data was filtered and the top {include_top_n_scoring}
      enriched regions per sample are retained.",
      " " = " "
    ))
  } else {
    # show error message independent of parameter show_messages
    options("rlib_message_verbosity" = "default")

    cli::cli_abort(c(
      "x" = "Given argument {.arg include_top_n_scoring} is not allowed.",
      "!" = "argument {.arg include_top_n_scoring} is
        '{.par {include_top_n_scoring}}'.",
      "i" = "Allowed values are NULL or single numeric value greater 1."
    ))
  }

  ### -----------------------------------------------------------------------###
  ### Return data
  ### -----------------------------------------------------------------------###

  return(data)
}
