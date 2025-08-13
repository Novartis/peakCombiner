#' Define which chromosomes are included
#'
#' @description
#' Retains only chromosomes that are in the provided vector. By not including
#' mitochondrial, sex, or non-classical chromosomes, genomic regions found on
#' these chromosomes can be removed. If set to 'NULL' (default), this step will
#' be skipped.
#'
#' @inheritParams filterRegions
#'
#' @return Data frame filtered by chromosome names based on the provided
#' parameters
#'
#' @noRd
#'
filterByChromosomeNames <- function(data,
                                    includeByChromosomeName = NULL,
                                    showMessages = TRUE) {
  
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
  ##
  ##
  ## check if input vector is numeric and if so change to character
  if (!is.character(includeByChromosomeName) &&
    !is.null(includeByChromosomeName)) {
    cli::cli_inform(c(
      "!" = "{.arg includeByChromosomeName} has the wrong class.",
      ">" = "{.arg includeByChromosomeName} is converted to
      {.cls character}."
    ))

    includeByChromosomeName <- as.character(includeByChromosomeName)
  }

  if (!is.null(includeByChromosomeName) &&
    !is.vector(includeByChromosomeName)) {
    # show error message independent of parameter show_messages
    options("rlib_message_verbosity" = "default")

    cli::cli_abort(c(
      "x" = "{.arg includeByChromosomeName} has to be a
        {.cls character} vector or {.val NULL}.",
      "!" = "Provided dataset is a
        {.cls {class(includeByChromosomeName)}}."
    ))
  } else if (length(includeByChromosomeName) > 1 &&
    any(is.null(includeByChromosomeName))) {
    # show error message independent of parameter show_messages
    options("rlib_message_verbosity" = "default")

    cli::cli_abort(c(
      "x" = "{.arg includeByChromosomeName} has a length of
      {length(includeByChromosomeName)} and contains {.val NULL}.",
      "i" = " Allowed values for {.arg includeByChromosomeName} are
      either a single {.val NULL} (length = 1) or vector with chromosome names
      (length >= 1)."
    ))
  } else if (any(is.na(includeByChromosomeName))) {
    # show error message independent of parameter show_messages
    options("rlib_message_verbosity" = "default")

    cli::cli_abort(c(
      "x" = "{.arg includeByChromosomeName} has a length of
      {length(includeByChromosomeName)} and contains {.val NA}.",
      "i" = " Allowed values for {.arg includeByChromosomeName} are
      either a single {.val NA} (length = 1) or vector with chromosome names
      (length >= 1)."
    ))
  } else if (is.null(includeByChromosomeName)) {
    cli::cli_inform(c(
      "i" = "The argument {.arg includeByChromosomeName} is {.val NULL}.",
      "v" = "No filtering for chromosome names in {.field chrom} is done.",
      " " = " "
    ))

    return(data)
  } else if (is.vector(includeByChromosomeName)) {
    # looks good so let's move on and check the goodness of the input vector

    cli::cli_inform(c(
      ">" = "The argument {.arg includeByChromosomeName} is a class
      {.cls character} of length {length(includeByChromosomeName)} and will
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
      includeByChromosomeName[!unique(includeByChromosomeName) %in%
        data_chr] |>
      unique()

    if (length(not_found_chr_names) > 0) {
      cli::cli_inform(c(
        "!" = "Input for {.arg includeByChromosomeName} contains values not
            found in the input data.",
        "i" = "The following chromosome name{?s} you entered {?is/are}
               not used: {.val {not_found_chr_names}}.",
        ">" = "{.emph Is {.val includeByChromosomeName} correctly defined?}"
      ))
    }
    rm(not_found_chr_names)

    # Get the chromosomes that were will be filtered out to report
    not_retained_chr_names <-
      data_chr[!data_chr %in% includeByChromosomeName] |> unique()

    retained_chr_names <-
      data_chr[data_chr %in% includeByChromosomeName] |> unique()

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
  
  ### -----------------------------------------------------------------------###
  ### Set message display back to default
  ### -----------------------------------------------------------------------###
  
  if (isFALSE(showMessages)) {
    options("rlib_message_verbosity" = "default")
  }

  return(data)
}

################################################################################
################################################################################
################################################################################

#' Define if and which blacklisted regions are excluded
#'
#' @description
#' A data frame or tibble can be provided listing the genomic regions to remove
#' (having `chrom`, `start`, and `end`  column names). If set to 'NULL'
#' (default), this step will be skipped.
#'
#' @inheritParams filterRegions
#'
#' @return Data frame filtered by blacklist based on the provided parameters.
#'
#' @noRd
#'
filterByBlacklist <- function(data,
                              excludeByBlacklist = NULL,
                              showMessages = TRUE) {
  
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
  ##
  

  ### -----------------------------------------------------------------------###
  ### Pre-Check up
  ### -----------------------------------------------------------------------###

  if (!exists("excludeByBlacklist")) {
    # show error message independent of parameter show_messages
    options("rlib_message_verbosity" = "default")

    cli::cli_abort(c(
      "x" = "Parameter {.arg excludeByBlacklist} doesn't exists.",
      "i" = "Allowed values are {.val {c('NULL')}} or
      a data frame with genomic coordinates in columns named {.field {c('chrom',
      'start', 'end')}}."
    ))
  }

  ### -----------------------------------------------------------------------###
  if (is.null(excludeByBlacklist)) {
    cli::cli_inform(c(
      "i" = "The argument {.arg excludeByBlacklist} is {.val NULL}.",
      "v" = "No filtering by blacklisted regions is done.",
      " " = " "
    ))

    return(data)
  } else if (is.data.frame(excludeByBlacklist)) {
    required_colnames_blacklist <- c("chrom", "start", "end")
    
    ## Check for correct colnames
    colnames(excludeByBlacklist) <- tolower(colnames(excludeByBlacklist))

    ## Check if data frame has chrom, start and end
    if (!all(required_colnames_blacklist %in% names(excludeByBlacklist))) {
      # show error message independent of parameter show_messages
      options("rlib_message_verbosity" = "default")

      cli::cli_abort(c(
        "x" = "{.arg excludeByBlacklist} is a data frame and misses
          some required column names.",
        "i" = "Required column names are:
        {.field {required_colnames_blacklist}}."
      ))
    }
    ## Check if chrom is character, start and end are numeric
    if (!is.character(excludeByBlacklist$chrom) &&
      !is.numeric(excludeByBlacklist$start) &&
      !is.numeric(excludeByBlacklist$end)) {
      # show error message independent of parameter show_messages
      options("rlib_message_verbosity" = "default")

      cli::cli_abort(c(
        "x" = "{.arg excludeByBlacklist} is a data frame and some columns
          are the wrong data type",
        "!" = "Required types for columns are: {.field chrom} as
        {.cls character}, {.field start} as {.cls numeric}, {.field end} as
        {.cls numeric} ."
      ))
    }

    cli::cli_inform(c(
      ">" = "User provied dataframe will be used for blacklist filtering."
    ))

    blacklist_data <- excludeByBlacklist
    
  } else if (inherits(excludeByBlacklist, "GRanges")) {
    cli::cli_inform(c(
      ">" = "User provied GenomicRanges object will be used for 
      blacklist filtering."
    ))
    
    blacklist_data <- dplyr::as_tibble(excludeByBlacklist) |>
      dplyr::rename(chrom = "seqnames") |>
      dplyr::select(-"width") |>
      dplyr::mutate(
        strand = ".",
        chrom = as.character(.data$chrom),
        strand = as.character(.data$strand)
      ) |>
      dplyr::ungroup()
    
    } else {
      cli::cli_abort(c(
        "x" = "The argument {.arg excludeByBlacklist} contains an 
        invalide value.",
        "i" = "The provided value is {.val {excludeByBlacklist}}"
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
  
  if (inherits(excludeByBlacklist, "GRanges")) {
    cli::cli_inform(c(
      "i" = "The argument {.arg excludeByBlacklist} is a class {.cls GRanges}.",
      ">" = "Using GenmoicRanges option for filtering."
    ))
    
    ## Do filtering: match of CHR between input and blacklist
    data <-
      data |>
      GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE) |>
      IRanges::subsetByOverlaps(
        excludeByBlacklist,
        invert = TRUE
      ) |> 
      suppressWarnings() |> #Recently added to solve warning
      tibble::as_tibble() |>
      dplyr::rename(chrom = "seqnames") |>
      dplyr::select(-"width") |>
      dplyr::mutate(
        strand = ".",
        chrom = as.character(.data$chrom),
        strand = as.character(.data$strand)
      ) |>
      dplyr::ungroup()
    
  } else {
    
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
        strand = ".",
        chrom = as.character(.data$chrom),
        strand = as.character(.data$strand)
      ) |>
      dplyr::ungroup()
  }
  
  cli::cli_inform(c(
    "v" = "Input data was filtered by blacklist.",
    " " = " "
  ))
  ### -----------------------------------------------------------------------###
  ### Set message display back to default
  ### -----------------------------------------------------------------------###
  
  if (isFALSE(showMessages)) {
    options("rlib_message_verbosity" = "default")
  }

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
#' `includeAboveScoreCutoff` parameter. If set to 'NULL' (default), this
#' step will be skipped.
#'
#' @inheritParams filterRegions
#'
#' @noRd
#'
filterBySignificance <- function(data,
                                 includeAboveScoreCutoff = NULL,
                                 showMessages = TRUE) {
  
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
  
  ##
  if (is.null(includeAboveScoreCutoff)) {
    cli::cli_inform(c(
      "i" = "The argument {.arg includeAboveScoreCutoff} is {.val NULL}.",
      "v" = "No filtering by {.field score} threshold is done.",
      " " = " "
    ))

    return(data)
  } else if (is.numeric(includeAboveScoreCutoff) &&
    length(includeAboveScoreCutoff) == 1) {
    cli::cli_inform(c(
      ">" = "Significance in {.field score} is filtered and all regions above
      {.val {includeAboveScoreCutoff}} will be retained."
    ))

    ## Format looks good, filter by the includeAboveScoreCutoff value

    input_rows <- nrow(data)

    data <-
      data |>
      dplyr::arrange(.data$sample_name, .data$score) |>
      dplyr::filter(.data$score >= includeAboveScoreCutoff) |>
      dplyr::ungroup()

    cli::cli_inform(c(
      "i" = "A total of {nrow(data)} of {input_rows} input regions are
      retained with value in {.field score} a above
      {includeAboveScoreCutoff}. ",
      "v" = "Input data was filtered to retain regions with a {.field score}
      above the defined threshold.",
      " " = " "
    ))
    rm(input_rows)
  } else {
    # show error message independent of parameter show_messages
    options("rlib_message_verbosity" = "default")

    cli::cli_abort(c(
      "x" = "Argument 'includeAboveScoreCutoff' has to be a
        {.cls numeric}.",
      "!" = "Provided argument is a
        {.cls {class(includeAboveScoreCutoff)}}"
    ))
  }
  
  ### -----------------------------------------------------------------------###
  ### Set message display back to default
  ### -----------------------------------------------------------------------###
  
  if (isFALSE(showMessages)) {
    options("rlib_message_verbosity" = "default")
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
#' discarded. Importantly, applying this filter retains `includeTopNScoring`
#' regions per sample, which means that the minimum enrichment levels may vary
#' between samples. Note that if multiple genomic regions have the same `score`
#' cutoff value, then all of those genomic regions are included. In this case,
#' the number of resulting regions retained may be a bit higher than the input
#' parameter. If set to 'NULL' (default), this step will be skipped.
#'
#' @inheritParams filterRegions
#'
#' @return Data frame filtered by top enriched regions based on the
#' provided parameters.
#'
#' @noRd
#'
filterByTopEnriched <- 
  function(data,
           includeTopNScoring = includeTopNScoring,
           showMessages = TRUE) {
  
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
  
    ##
  ##
  if (is.null(includeTopNScoring)) {
    cli::cli_inform(c(
      "i" = "The argument {.arg includeTopNScoring} is {.val NULL}.",
      "v" = "No top enriched regions were selected. All input regions are
      retained.",
      " " = " "
    ))

    return(data)
  } else if (is.numeric(includeTopNScoring) &&
    includeTopNScoring > 0) {
    ### ---------------------------------------------------------------------###

    cli::cli_inform(c(
      "i" = "The argument {.arg includeTopNScoring} extracted the the top
      {.num {includeTopNScoring}} regions by {.field score} per sample (based
      on the values in {.field sample_name}).",
      ">" = "The top enriched {.num {includeTopNScoring}} regions per sample
      will be retained."
    ))

    too_few_regions_left <-
      data |>
      dplyr::group_by(.data$sample_name) |>
      dplyr::summarise(counts = dplyr::n(), .groups = "drop") |>
      dplyr::filter(.data$counts < includeTopNScoring) |>
      dplyr::pull(.data$sample_name)

    if (length(too_few_regions_left) > 0) {
      cli::cli_inform(c(
        "i" = "The argument {.arg includeTopNScoring} was defined as
        {includeTopNScoring}.",
        ">" = "The following {.val sample_names} contain less regions then
        defined by {.arg includeTopNScoring}: {too_few_regions_left}",
        "!" = "No genomic regions will be removed for such samples."
      ))
    }

    ### ---------------------------------------------------------------------###

    data <-
      data |>
      dplyr::group_by(.data$sample_name) |>
      dplyr::top_n(n = !!includeTopNScoring, wt = .data$score) |>
      dplyr::ungroup()

    cli::cli_inform(c(
      "v" = "Input data was filtered and the top {includeTopNScoring}
      enriched regions per sample are retained.",
      " " = " "
    ))
  } else {
    # show error message independent of parameter show_messages
    options("rlib_message_verbosity" = "default")

    cli::cli_abort(c(
      "x" = "Given argument {.arg includeTopNScoring} is not allowed.",
      "!" = "argument {.arg includeTopNScoring} is
        '{.par {includeTopNScoring}}'.",
      "i" = "Allowed values are NULL or single numeric value greater 1."
    ))
  }
    
  ### -----------------------------------------------------------------------###
  ### Set message display back to default
  ### -----------------------------------------------------------------------###

  if (isFALSE(showMessages)) {
    options("rlib_message_verbosity" = "default")
  }

  ### -----------------------------------------------------------------------###
  ### Return data
  ### -----------------------------------------------------------------------###

  return(data)
}
