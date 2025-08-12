#' Modifies genomic regions by centering and then expanding them
#'
#' @description
#' [peakCombiner::centerExpandRegions()] is an optional step that re-defines the
#' genomic regions by expanding them from their center. The center information
#' has to be stored in the input data column `center`, while the information for
#' the expansion can either be user provided or input data derived. The accepted
#' input is a data frame created from [peakCombiner::prepareInputRegions()].
#' Please see [peakCombiner::prepareInputRegions()] for more details.
#'
#' @details
#' This is an optional function that resizes the genomic regions based on the
#' input peakCombiner standard data frame and the options you select. An
#' expected input data foam contains the following columns with the names:
#' `chrom`, `start`, `end`, `name`, `score`, `strand`, `center`, `sample_name`.
#' Such a data frame is created by the script
#' [peakCombiner::prepareInputRegions]. This step is useful if you want all of
#' your peaks to be the same size for your downstream analyses. In addition, if
#' you want to use the "summit" information, normally obtained by some peak
#' callers (e.g., Macs2), this function allows you to automatically center your
#' regions of interest on these summits. This enables you to capture
#' information about the most important region within a genomic region (e.g.,
#' TF-binding site or highest peak) and put that region in the center of your
#' downstream analyses (e.g., applicable to motif-finding or "heatmaps"
#' summarizing multiple genomic regions).
#'
#' There are two concepts that are relevant for
#' [peakCombiner::centerExpandRegions]: how to define the center, and how much
#' to expand from the center.
#'
#' ## How to define the center?
#'
#' When you prepared your input regions, it is recommended to use the function
#' [peakCombiner::prepareInputRegions] provided by this package. This pre-
#' populated the `center` column with the absolute genomic coordinate of the
#' center of the peak region. You can either choose to define the center by
#' using pre-defined summit information (e.g., obtained from a peak caller like
#' MACS2) or re-compute the arithmetic mean and save that value in the column
#' `center`. (For details see the help for
#' [peakCombiner::prepareInputRegions()]).
#'
#' ## How much to expand from the center
#' You can choose to expand the genomic region from the center either
#' symmetrically or asymmetrically (different lengths before and after the
#' center position).
#'
#' In the symmetrical case, if you want to choose the size of your genomic
#' region based on the input data, this function can also calculate the median
#' peak size across all of your genomic regions and use that value (`expandBy`
#' = NULL). Alternatively, the user is free to provide a numeric vector to
#' define the expansion. A numeric vector with one value is used to
#' symmetrically expand, while a vector with two values allows to expand
#' asymmetrically.
#'
#'
#' @param data        PeakCombiner data frame structure with required columns
#'                      named `chrom`, `start`, `end`, `name`,
#'                      `score`, `strand`, `center`, `sample_name`. Additional
#'                      columns will be maintained.
#' @param centerBy   Allowed values are 'center_column' (default) or
#'                    'midpoint'.
#' * 'center_column' uses the value stored in the column `center` to center.
#' * 'midpoint' replaces the value stored in the column `center` based on the 
#' [GenomicRanges::resize()] followed by the expansion from based on the user 
#' input using [GenomicRanges::promoters()] to allow symmetric and asymmetic 
#' expansion. Note that strand information, if provided is maintained.
#'
#' @param expandBy   Allowed values a numeric vector of length 1 or 2,
#'                      or 'NULL' (default).
#' * The value from the numeric vector of length 1
#'                          is expanded in both directions from center to define
#'                          the genomic region.
#'                          Thus, the size of the resulting genomic region is 2x
#'                          the provided value + 1 (for the center coordinate).
#' * The value of the numeric vector of length 2
#'                          subtracts the first value from the center and adds
#'                          the second value to the center to define the genomic
#'                          region. Thus, the size of the genomic regions is
#'                          the sum of the first value + the second value
#'                          + 1 (for the center coordinate).
#' * 'NULL' allows for data-driven definition of the
#'                          `expandBy` value. It calculates the median
#'                          genomic region size of the input data and uses this
#'                          value like a length 1 numeric vector for expansion.
#' @param genome      Character value to define the matching genome reference to 
#'                      the input data. Default value is NA. Allows values are 
#'                      based on GenomicRanges supported genomes like "GRCh38", 
#'                      "GRCh38.p13", "Amel_HAv3.1", "WBcel235", "TAIR10.1", 
#'                      "hg38", "mm10", "rn6", "bosTau9", "canFam3", "musFur1", 
#'                      "galGal6","dm6", "ce11", and "sacCer3". Please see also
#'                      help for [GenomeInfoDb::Seqinfo()] for more details.
#'
#' @param trim_start  Logical value of TRUE or FALSE (default). If TRUE, and 
#'                      no valid reference genome are provided in `genome`, 
#'                      resulting genomic results with negative starting 
#'                      coordinates will be set to 1.
#'            
#' @param outputFormat Character value to define format of output object. 
#'                      Accepted values are "GenomicRanges" (default), "tibble" 
#'                      or "data.frame".  
#'
#' @param showMessages Logical value of TRUE (default) or FALSE. Defines if
#'                      info messages are displayed or not.
#'
#'
#' @return A tibble with the columns `chrom`, `start`, `end`, `name`, `score`,
#' `strand`, `center`, `sample_name`. The definitions of these columns are
#' described in full in the [peakCombiner::prepareInputRegions] Details.
#' Use as input for functions [peakCombiner::filterRegions()] and
#' [peakCombiner::combineRegions()].
#'
#' @export
#'
#' @importFrom rlang .data
#' @import tidyr
#' @import here
#'
#' @examples
#' # Load in and prepare a an accepted tibble
#' utils::data(syn_data_bed)
#'
#' # Prepare input data
#' data_prepared <- prepareInputRegions(
#'   data = syn_data_bed,
#'   outputFormat = "tibble",
#'   showMessages = TRUE
#' )
#' # Run center and expand
#' data_center_expand <- centerExpandRegions(
#'   data = data_prepared,
#'   centerBy = "center_column",
#'   expandBy = NULL,
#'   outputFormat = "tibble",
#'   showMessages = TRUE
#' )
#'
#' data_center_expand
#'
#' # You can choose to use the midpoint and predefined values to expand
#'
#' data_center_expand <- centerExpandRegions(
#'   data = data_prepared,
#'   centerBy = "midpoint",
#'   expandBy = c(100, 600),
#'   outputFormat = "tibble",
#'   showMessages = FALSE
#' )
#'
#' data_center_expand
#'
centerExpandRegions <- function(data,
                                  centerBy = "center_column",
                                  expandBy = NULL,
                                  genome = NA,
                                  trim_start = TRUE,
                                  outputFormat = "GenomicRanges",
                                  showMessages = TRUE) {
  
  ### -----------------------------------------------------------------------###
  ### Allowed genomes from GenomicRanges
  ### -----------------------------------------------------------------------###

  gr_genome_seqinfo <- c("GRCh38", 
                         "GRCh38.p13", 
                         "Amel_HAv3.1", 
                         "WBcel235", 
                         "TAIR10.1",
                         "hg38", 
                         "mm10", 
                         "rn6", 
                         "bosTau9", 
                         "canFam3", 
                         "musFur1", 
                         "galGal6",
                         "dm6",  
                         "ce11", 
                         "sacCer3")
  
  column_order <- colnames(data)
  
  ### -----------------------------------------------------------------------###
  ### Prepare parameters
  ### -----------------------------------------------------------------------###
  center_values <- c("center_column", "midpoint")
  
  ## Check parameter value correctness and calculate if needed
  expansion_value <- defineExpansion(
    data = data,
    expandBy = expandBy
  )
  
  ## Calculate the values to expand the regions
  length_expansion_value <- length(expansion_value)
  expand_1 <- expansion_value[1]
  expand_2 <- expansion_value[length_expansion_value]
  
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
    
  } else {
    cli::cli_inform(c(
      "i" = "Input data {.arg data} has no assigned genome."
    ))  
    input_file_genome <- NA
  }
  
  
  if (!is.na(input_file_genome) & is.na(genome)) {
    cli::cli_inform(c(
    "i" = "Input data {.arg data} assigned genome ({.val {input_file_genome}})
    is used for trimming."
    ))
    
    genome_used <- input_file_genome
      
    } else if (!is.na(input_file_genome) & genome %in% gr_genome_seqinfo){
      cli::cli_inform(c(
        "!" = "Input data {.arg data} assigned genome 
      ({.val {input_file_genome}}) is replaced by provided {.arg genome}.",
        "i" = "Argument {.arg genome} set to {.val {genome}} is used."
      ))
      
      genome_used <- genome
      
    } else if (is.na(input_file_genome) & genome %in% gr_genome_seqinfo) {
      cli::cli_inform(c(
        "!" = "Input data {.arg data} has no assigned genome 
      ({.val {input_file_genome}}).",
        "i" = "Argument {.arg genome} set to {.val {genome}} is used."
      ))
      
      genome_used <- genome
      
    } else if (is.na(input_file_genome) & is.na(genome)) {
      cli::cli_inform(c(
        "!" = "Input data {.arg data} has no assigned genome 
      ({.val {input_file_genome}}).",
        "i" = "Argument {.arg genome} set to {.val {genome}}.",
        "i" = "Only start will be trimmed."
      ))
      
      genome_used <- NA
      
    } else {
      cli::cli_abort(c(
        "x" = "Argument {.arg genome} has to be one of the following
      values: {.val {gr_genome_seqinfo}}.",
        "i" = "Provided value is {.val {genome}}."
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
      "!" = "Provided input {.arg data} is a class {.cls GRanges}.",
      ">" = "Start converting and preparing data."
    ))
    
    input_seqinfo <- GenomeInfoDb::seqinfo(data)
    
    data <-
      tibble::as_tibble(data) |>
      dplyr::rename(chrom = "seqnames") |>
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
    
  } else {
    # show error independend of showMessages
    options("rlib_message_verbosity" = "default")
    
    cli::cli_abort(c(
      "x" = "Provide input {.arg data} does not have the required format.",
      "!" = "Please check your column names in {.arg data}."
    ))
  }
  
  ### -----------------------------------------------------------------------###
  ### Trim start coordinates only 
  ### -----------------------------------------------------------------------###
  
  if (!is.logical(trim_start)) {
    cli::cli_abort(c(
      "x" = "Argument {.arg trim_start} has to be {.cls logical}."
    ))
  } else {
    # show error message independent of parameter trim_start
    options("rlib_message_verbosity" = "default")
    
    cli::cli_inform(c(
      "i" = "Argument {.arg trim_start} is {.val {trim_start}}."
    ))
  }
 
  ### -----------------------------------------------------------------------###
  ### Check input parameters
  ### -----------------------------------------------------------------------###

  if (is.null(centerBy)) {
    # show error message independent of parameter showMessages
    options("rlib_message_verbosity" = "default")

    cli::cli_abort(c(
      "x" = "{.arg centerBy} has to be {.val center_column} or
      {.val midpoint}.",
      "i" = "{.arg centerBy} is {.val NULL}."
    ))
  } else if (length(centerBy) != 1) {
    # show error message independent of parameter showMessages
    options("rlib_message_verbosity" = "default")

    cli::cli_abort(c(
      "x" = "{.arg centerBy} has a length of {length(centerBy)}.",
      "i" = "{.arg centerBy} allowed length is 1."
    ))
  } else if (!tolower(centerBy) %in% center_values) {
    # show error message independent of parameter showMessages
    options("rlib_message_verbosity" = "default")

    cli::cli_abort(c(
      "x" = "{.arg centerBy} has to be {.val center_column} or
      {.val midpoint}.",
      "i" = "{.arg centerBy} is {.val {centerBy}}."
    ))
  } else if (tolower(centerBy) %in% center_values) {
    ## good values!
    centerBy <- tolower(centerBy)
  } else {
    # show error message independent of parameter showMessages
    options("rlib_message_verbosity" = "default")

    cli::cli_abort(c(
      "x" = "{.arg centerBy} has to be {.val center_column} or
      {.val midpoint}.",
      "i" = "{.arg centerBy} is {.val NULL}."
    ))
  }

  ## Check the validity of the peakCombiner input data format

  data <- checkDataStructure(
    data = data
  )

  ### -----------------------------------------------------------------------###
  ### Center and expand
  ### -----------------------------------------------------------------------###
  
  if (centerBy == "center_column") {
    
    cli::cli_inform(c(
      ">" = "Starting with expanding genomic regions from the column 
      {.field center}."
    ))
    
  cli::cli_inform(c(
    ">" = "Genomic regions will be centered and expanded.",
    "i" = "Used {.field genome} for trimming is {.val {genome_used}}."
  ))
  
  if(!is.na(genome_used)) {
      # Use provided genome
      gr_genome <- GenomeInfoDb::Seqinfo(genome = genome_used)
      gr_seqlevels <- GenomeInfoDb::seqlevels(gr_genome)
      
      # Filter not matching chromsomes
      data_filtered <- data |>
        dplyr::filter(.data$chrom %in% gr_seqlevels)
      
      data_gr <- data_filtered |>
        dplyr::mutate(center = round(.data$center, 0),
                      start = .data$center,
                      end = .data$center + 1)
      
      # Convert to GRanges
      gr_resized <- data_gr |>
        GenomicRanges::makeGRangesFromDataFrame(
        keep.extra.columns = TRUE,
        seqinfo = gr_genome,
        ignore.strand = FALSE
      )

      data_center_expand <- GenomicRanges::trim(
        GenomicRanges::promoters(
          gr_resized, 
          upstream = expand_1, 
          downstream = expand_2,
          use.names = TRUE)) |>
        suppressWarnings()

    } else {
      # Convert to GRanges
      data_gr <- data |>
        dplyr::mutate(center = round(.data$center, 0),
                      start = .data$center,
                      end = .data$center) 
      
      gr_resized <- data_gr |>
        GenomicRanges:: makeGRangesFromDataFrame(
          keep.extra.columns = TRUE,
          seqinfo = NULL,
          ignore.strand = FALSE
        )
      
      data_center_expand <- GenomicRanges::promoters(gr_resized, 
                                           upstream = expand_1, 
                                           downstream = expand_2, 
                                           use.names = TRUE)
    }
    
    # Show information about expansion
    cli::cli_inform(c(
      ">" = "Expanding genomic regions from the column {.field center} by
      {.val {expand_1}} before and {.val {expand_2}} after the center."
    ))
    
    } else if (centerBy == "midpoint") {
    cli::cli_inform(c(
      ">" = "Starting with defining the {.field center} as the midpoint of the 
      regions from the {.field start} and {.field end} coordinates."
    ))
      
      if(!is.na(genome)) {
        cli::cli_inform(c(
          ">" = "Using {.field genome} {.val {genome}} to assign genome."
        ))
        
        gr_genome <- GenomeInfoDb::Seqinfo(genome = genome)
        gr_seqlevels <- GenomeInfoDb::seqlevels(gr_genome)
        
        # Filter not matching chromsomes
        data_filtered <- data |>
          dplyr::filter(.data$chrom %in% gr_seqlevels)
        
        # Convert to GRanges
        data_gr <- GenomicRanges::makeGRangesFromDataFrame(
          df = data_filtered,
          keep.extra.columns = TRUE,
          seqinfo = gr_genome,
          ignore.strand = FALSE
        )
        
        # Expand with reference genome and trimming
        gr_resized <- GenomicRanges::resize(data_gr, 
                                            width = 2, 
                                            fix = "center",
                                            ignore.strand=FALSE)
        
        data_center_expand <- GenomicRanges::trim(
          GenomicRanges::promoters(gr_resized, 
                                   upstream = expand_1, 
                                   downstream = expand_2, 
                                   use.names = TRUE)) |>
          suppressWarnings()
        
      } else {
        cli::cli_inform(c(
          ">" = "Option {.field genome} is empty, so no genome is assigned."
        ))
        
        # Convert to GRanges
        data_gr <- GenomicRanges::makeGRangesFromDataFrame(
          df = data,
          keep.extra.columns = TRUE,
          seqinfo = NULL,
          ignore.strand = FALSE
        )
        # Resize to center
        gr_resized <- GenomicRanges::resize(data_gr, width = 2, fix = "center")
        
        # Expand without reference genome, no trimming
        data_center_expand <- GenomicRanges::promoters(gr_resized, 
                                                       upstream = expand_1,
                                                       downstream = expand_2,
                                                       use.names = TRUE)
        
      }
    }
  
  # Convert back to tibble
  data_center_expand <- dplyr::as_tibble(data_center_expand) |>
    dplyr::rename(chrom = "seqnames" ) |>
    dplyr::select(column_order[c(1:3, 6, 4:5, 7:8)]) |>
    dplyr::mutate(strand = ".",
                  chrom = as.character(.data$chrom)) |>
    dplyr::select(any_of(column_order))
  
  
  # View result
  print(data_center_expand)
  
  cli::cli_inform(c(
    "v" = "Genomic regions were successfully centered and expanded.",
    " " = " "
  ))

  
  if(isTRUE(trim_start) & is.na(genome)) {
    cli::cli_inform(c(
      "i" = "Argument {.arg trim_start} is set to {.val TRUE}.",
      "i" = "Atgument {.arg genome} is set to NA.",
      ">" = "Trimming start coordinates of resulting genomic regions."
    ))
   
     if (any(data_center_expand$start < 1)) {
      neg_starts <-
        data_center_expand |>
        dplyr::filter(.data$start < 1) |>
        dplyr::pull(.data$name)
      
      cli::cli_inform(c(
        "i" = "Some newly-defined genomic regions have a {.field start}
      coordinate below {.val 1}.",
        ">" = "Values of {.field name} for these  site{?s}: {.val {neg_starts}}."
      ))
      
      data_center_expand <-
        data_center_expand |>
        dplyr::mutate(start = ifelse((.data$start < 1), 1, .data$start))
      
      cli::cli_inform(c(
        "v" = "These genomic regions were trimmed to get {.field start}
      coordinate {.val 1}."
      ))
      
      rm(neg_starts)
    }
  }

  ### -----------------------------------------------------------------------###
  ### Return data frame
  ### -----------------------------------------------------------------------###

  cli::cli_inform(c(
    "v" = "Genomic regions were successfully centered and expanded.",
    " " = " "
  ))
  
  ### -----------------------------------------------------------------------###
  ### Adjust output format
  ### -----------------------------------------------------------------------###
  
  if (outputFormat %in% c("GenomicRanges", "GRanges")) {
    #if(exists("input_seqinfo")) {
    if(!is.na(genome_used)) {
    cli::cli_inform(c(
        "i" = "Output format is set to {.val {outputFormat}}.",
        "i" = "Assigning {.val {genome_used}} genome annotation to output. ")
      )
      
      data_center_expand <- 
        data_center_expand |>
        GenomicRanges::makeGRangesFromDataFrame(
          keep.extra.columns = TRUE,
          seqinfo = gr_genome, 
          starts.in.df.are.0based = FALSE
        )
    } else {
      cli::cli_inform(c(
        "i" = "Output format is set to {.val {outputFormat}}.",
        "i" = "No input genome annotation assigned to ouutput. ")
      )
      data_center_expand <- 
        data_center_expand |>
        GenomicRanges::makeGRangesFromDataFrame(
          keep.extra.columns = TRUE,
          starts.in.df.are.0based = FALSE
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

  return(data_center_expand)
}
