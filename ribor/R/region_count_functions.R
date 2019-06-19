#' Retrieves the region counts from a .ribo file
#'
#' \code{\link{get_region_counts}} will return the particular region counts
#' of any subset of regions for a given set of experiments.
#'
#' This function will return a data.table of the counts at each specified region
#' for each specified experiment. The region options are "UTR5", "UTR5J", "CDS",
#' "UTR3J", and "UTR3". The user can specify any subset of regions in the form of a vector,
#' a list, or a single string if only one region is required.
#'
#' The dimensions of the returned data table depend on the parameters
#' range.lower, range.upper, length, and transcript.
#'
#' The param 'length' condenses the read lengths together.
#' When length is true and transcript is false, the
#' data table presents information for each transcript across
#' all of the read lengths. That is, each transcript has a value
#' that is the sum of all of the counts across every read length.
#' As a result, information about the transcript at each specific
#' read length is lost.
#'
#' The param 'transcript' condenses the transcripts together.
#' When transcript is true and length.lengths is false, the data
#' table presents information at each read length between range.lower and
#' range.upper inclusive. That is, each separate read length denotes the
#' sum of counts from every transcript. As a result, information about the
#' counts of each individual transcript is lost.
#'
#' If both 'length' and 'transcript' are true, then the resulting
#' data table prints out one row for each experiment. This provides the metagene
#' information across all transcripts and all reads in a given experiment. If both
#' 'length' and 'transcript' are true, no calculations are done to the data,
#' all information is preserved for both the read length and the transcript.
#' The data table would just present the entire stored raw data
#' from the read length 'range.lower' to the read length 'range.upper' which in most
#' cases would result in a slow run time with a massive data.table returned.
#'
#' @examples
#' #generate the ribo object
#' file.path <- system.file("extdata", "sample.ribo", package = "ribor")
#' sample <- ribo(file.path)
#'
#' #specify the regions and experiments of interest
#' regions <- c("UTR5", "UTR5J", "CDS", "UTR3J", "UTR3")
#' experiments <- c("Hela_1", "Hela_2", "WT_1")
#'
#' #obtains the region counts at each individual read length, summed across every transcript
#' region.counts <- get_region_counts(sample,
#'                                    regions,
#'                                    range.lower = 2,
#'                                    range.upper = 5,
#'                                    length = FALSE,
#'                                    transcript = TRUE,
#'                                    experiments = )
#'
#' @param ribo.object S3 "ribo" class object
#' @param regions Specific region of interest
#' @param range.lower Lower bound of the read length
#' @param range.upper Upper bound of the read length
#' @param length Option to condense the read lengths together, preserve the transcripts
#' @param transcript Option to condense the transcripts together, preserve the read lengths
#' @param experiments List of experiment names
#' @return A data table of the region counts
#' @importFrom rhdf5 h5read
#' @importFrom data.table data.table
#' @export
get_region_counts <- function(ribo.object,
                              regions,
                              range.lower,
                              range.upper,
                              length = TRUE,
                              transcript = TRUE,
                              experiments = get_experiments(ribo.object)) {
  if (typeof(regions) == 'character') {
    regions <- toupper(as.vector(regions))
  }

  range.info <- c(range.lower = range.lower,
                  range.upper = range.upper)

  check_rc_input(ribo.object,
                 regions,
                 range.info,
                 experiments)

  conditions <- c(transcript = transcript,
                     length = length)

  ribo.experiments <- get_experiments(ribo.object)
  handle <- ribo.object$handle
  matched.experiments <- intersect(experiments, ribo.experiments)
  total.experiments <- length(matched.experiments)

  region.list <- c(1, 2, 3, 4, 5)
  names(region.list) <- c("UTR5", "UTR5J", "CDS", "UTR3J", "UTR3")

  ref.names    <- get_reference_names(ribo.object)
  ref.length   <- length(get_reference_names(ribo.object))
  read.lengths <- get_read_lengths(ribo.object)
  ribo.min     <- read.lengths[1]
  read.length.range <- range.upper - range.lower + 1

  result <- data.table()

  for (region in regions) {
    column.index <- region.list[[toupper(region)]]
    region.data <- determine_matrix_size(conditions,
                                         ref.length,
                                         1,
                                         total.experiments,
                                         read.length.range)

    colnames(region.data) <- "region.count"

    region.info <- list(ribo.min            = ribo.min,
                        range.info          = range.info,
                        ref.names           = ref.names,
                        total.experiments   = total.experiments,
                        matched.experiments = matched.experiments,
                        ref.length          = ref.length,
                        column.index        = column.index,
                        conditions          = conditions)

    current.data <- fill_region(handle, region.info, region.data)
    result       <- rbind(result, current.data)
  }

  #calculates how many rows that each region has
  rep <- nrow(result)/length(regions)

  return(cbind(result, region = rep(unlist(regions), each = rep)))
}


fill_region <- function (handle,
                         region.info,
                         region.data) {
  #helper function that fills information for a single regions
  ribo.min            <- region.info[["ribo.min"]]
  range.lower         <- region.info[["range.info"]][["range.lower"]]
  range.upper         <- region.info[["range.info"]][["range.upper"]]
  ref.names           <- region.info[["ref.names"]]
  total.experiments   <- region.info[["total.experiments"]]
  matched.experiments <- region.info[["matched.experiments"]]
  ref.length          <- region.info[["ref.length"]]
  column.index        <- region.info[["column.index"]]
  transcript          <- region.info[["conditions"]][["transcript"]]
  length              <- region.info[["conditions"]][["length"]]

  for (i in 1:total.experiments) {
    experiment <- matched.experiments[i]
    path <- paste("/experiments/", experiment,
                  "/region_counts/region_counts", sep ="")
    for (current.length in (range.lower:range.upper)) {
      offset <- ref.length * (current.length - ribo.min) + 1
      row.index    <- c(offset: (offset + ref.length - 1))
      region.count <- t(h5read(handle, path,
                               index = list(column.index, row.index)))
      fill.params <- list(index = i,
                          conditions = region.info[["conditions"]],
                          ref.length = ref.length,
                          result = region.data,
                          current.length = current.length,
                          range.info = region.info[["range.info"]],
                          data = region.count)

      region.data <- fill_matrix(fill.params)
    }
  }
  return(make_datatable(ref.names,
                        transcript,
                        length,
                        matched.experiments,
                        range.lower,
                        range.upper,
                        region.data))
}


check_rc_input <- function(ribo.object,
                           regions,
                           range.info,
                           experiments) {
  #helper function that checks for valid parameters given by the user
  #calls error messages on any incorrect parameters
  region.options <- c("UTR5", "UTR5J", "CDS", "UTR3J", "UTR3")

  if (typeof(regions) != 'list' & typeof(regions) != "character") {
    stop("Please specify the regions as a single string, a vector, or a list.")
  }

  for (region in regions) {
    if (!region %in% region.options) {
      stop("'", region, "'", " is not an option. Please indicate the region with
         one of the following: 'UTR5', 'UTR5J, 'CDS', 'UTR3J', 'UTR3'")
    }
  }

  range.lower <- range.info[["range.lower"]]
  range.upper <- range.info[["range.upper"]]
  check_lengths(ribo.object, range.lower, range.upper)
  check_experiments(ribo.object, experiments)
}

#' Plots the length distribution
#'
#' The function \code{\link{plot_region_counts}} can take either a data.table
#' or a "ribo" object to generate a line graph of the length distributions from
#' range.lower to range.upper.
#'
#' The param 'percentage' will plot the percentages of each length relative
#' to the total sum of the read length range provided by param 'range.lower'
#' and 'range.upper'. When percentage is set to FALSE, the total count of each
#' read length is plotted.
#'
#' When given a "ribo" object, \code{\link{plot_length_distribution}} calls
#' \code{\link{get_region_counts}} to retrieve the necessary information
#' for plotting. This option is in the case that a data.table of the
#' region count information is not required.
#'
#' The user can instead provide a data.table with the same structure as the
#' output of the \code{\link{get_region_counts}} function where the 'transcript'
#' parameter is set to FALSE and 'length' parameters is the default value of
#' TRUE. This also means that the remaining parameters of the
#' \code{\link{plot_length_distribution}} function are not necessary. The run
#' time becomes substantially faster when \code{\link{plot_region_counts}} is
#' given the direct data.table to plot. Note that there is no manipulation by
#' this function on the data.table, making this input option more error prone.
#'
#' @examples
#' #ribo object use case
#'
#' #generate the ribo object
#' file.path <- system.file("extdata", "sample.ribo", package = "ribor")
#' sample <- ribo(file.path)
#'
#' #specify the regions and experiments of interest
#' regions <- c("UTR5", "UTR5J", "CDS", "UTR3J", "UTR3")
#' experiments <- c("Hela_1", "Hela_2", "wt_1")
#'
#' plot_length_distribution(sample,
#'                          region = "CDS",
#'                          range.lower = 2,
#'                          range.upper = 5,
#'                          experiments = experiments,
#'                          percentage = TRUE)
#'
#'
#' #data.table use case
#' #obtains the region counts at each individual read length, summed across every transcript
#' region.counts <- get_region_counts(sample,
#'                                    regions,
#'                                    range.lower = 2,
#'                                    range.upper = 5,
#'                                    length = FALSE,
#'                                    transcript = TRUE)
#'
#' #the param 'length' must be set to FALSE and param 'transcript' must be set
#' #to TRUE to use a data.table
#' plot_length_distribution(region.counts)
#'
#'
#' @seealso \code{\link{get_region_counts}} to generate a data.table that can
#' be provided as input,
#' \code{\link{ribo}} to create a ribo.object that can be provided as input
#' @param x either an S3 "ribo" object or a data.table
#' @param region the region of interest
#' @param range.lower a lower bounds for a read length range
#' @param range.upper an upper bounds for a read length range
#' @param experiments a list of experiment names
#' @param percentage logical value that, if TRUE, presents the data as a percentage
#' @param title a title for the generated plot
#' @importFrom data.table is.data.table
#' @importFrom tidyr gather
#' @importFrom dplyr select
#' @importFrom ggplot2 ggplot geom_line theme_bw theme labs
#' @importFrom rlang sym
#' @importFrom stats aggregate
#' @importFrom rlang .data
#' @export
plot_length_distribution<- function(x,
                                    region,
                                    range.lower,
                                    range.upper,
                                    experiments,
                                    percentage = TRUE,
                                    title = "Length Distribution"){
  is.ribo <- check_ribo(x, stop = FALSE)
  if(is.ribo) {
    missing.ranges <- missing(range.lower) || missing(range.upper)
    if (missing.ranges) {
      stop("Please indicate the 'range.lower' and 'range.upper' parameters.")
    } else if (missing(region)) {
      stop("Please indicate the 'region' parameters.")
    } else if (missing(experiments)) {
      experiments = get_experiments(x)
    } else if (missing(range.lower)) {
      range.lower <- get_read_lengths(x)[1]
    } else if (missing(range.upper)) {
      range.upper <- get_read_lengths(x)[2]
    }

    x <- get_region_counts(x,
                           region,
                           range.lower,
                           range.upper,
                           length = FALSE,
                           transcript = TRUE,
                           experiments = experiments)
  } else if (is.data.table(x)){
    col.names <- c("experiment", "length", "region.count", "region")
    mismatch <- !all(names(x) == col.names) || typeof(x[[1]]) != "character"||
    typeof(x[[2]]) != "integer"  || typeof(x[[3]]) != "double" ||
    typeof(x[[4]]) != "character" || ncol(x) != 4
    if (mismatch) {
      stop("Please make sure that the data table is of the correct format.")
    }
  } else {
    stop("Please make sure that param 'x' is either a data.table or a ribo object.")
  }

  y.axis <- "Count"
  y.value <- sym("region.count")

  if (percentage) {
    total <- aggregate(x$region.count, by = list(experiment= x$experiment), FUN=sum)
    x <- x %>%  left_join(total, by= "experiment") %>%
      mutate(fraction = .data$region.count/x) %>%
      select(-x)
    y.axis <- "Frequency"
    y.value <- sym("fraction")
  }

  ggplot(x,
         aes(x     = length,
             y     = !!y.value,
             color = .data$experiment)) +
    geom_line() +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5)) +
    labs(title = title, x = "Read Length", y = y.axis, color = "Experiment")
}

#' Plots the region counts of UTR5, CDS, and UTR3
#'
#' The function \code{\link{plot_region_counts}} can take either a data.table
#' or a "ribo" object to generate the a stacked bar plot of proportions that
#' correspond to the "UTR5", "CDS", and "UTR3" regions.
#'
#' When given a "ribo" object, \code{\link{plot_region_counts}} calls
#' \code{\link{get_region_counts}} to retrieve the necessary information
#' for plotting. This option is in the case that a data.table of the
#' region count information is not required.
#'
#' The user can instead provide a data.table with the same structure as the
#' output of the \code{\link{get_region_counts}} function where the 'transcript'
#' and 'length' parameters are the default values of TRUE. This also means that
#' the remaining parameters of the \code{\link{plot_region_counts}} function are not necessary.
#' The run time becomes substantially faster when \code{\link{plot_region_counts}} is given
#' the direct data.table to plot. Note that there is no manipulation by this function on the
#' data.table, making this input option more error prone.
#'
#' @examples
#' #ribo object use case
#' #generate the ribo object
#' file.path <- system.file("extdata", "sample.ribo", package = "ribor")
#' sample <- ribo(file.path)
#'
#' #specify the regions and experiments of interest
#' regions <- c("UTR5", "CDS", "UTR3")
#' experiments <- c("Hela_1", "Hela_2", "WT_1")
#'
#' plot_region_counts(sample,
#'                    range.lower = 2,
#'                    range.upper = 5,
#'                    experiments)
#'
#' #data.table use case
#' #obtains the region counts at each individual read length, summed across every transcript
#' region.counts <- get_region_counts(sample,
#'                                    regions,
#'                                    range.lower = 2,
#'                                    range.upper = 5,
#'                                    length = TRUE,
#'                                    transcript = TRUE)
#'
#' #the params 'length' and 'transcript' must be set to true to use a data.table
#' plot_region_counts(region.counts)
#'
#' @seealso \code{\link{get_region_counts}} to generate a data.table that can be provided as input,
#' \code{\link{ribo}} to create a ribo.object that can be provided as input
#'
#' @param x either an S3 "ribo" object or a data.table
#' @param range.lower a lower bounds for a read length range
#' @param range.upper an upper bounds for a read length range
#' @param experiments a list of experiment names
#' @param title a title for the generated plot
#' @importFrom data.table is.data.table
#' @importFrom dplyr left_join mutate %>% group_by summarize arrange desc
#' @importFrom tidyr gather
#' @importFrom rlang .data
#' @importFrom ggplot2 ggplot geom_col theme_bw theme ggtitle coord_flip theme
#' @importFrom ggplot2 labs scale_fill_discrete
#' @export
plot_region_counts <- function(x,
                               range.lower,
                               range.upper,
                               experiments,
                               title = "Region Counts") {
  is.ribo <- check_ribo(x, stop = FALSE)
  regions <- c("UTR5", "CDS", "UTR3")
  all.regions <- data.table()
  if (is.ribo) {
    if (missing(experiments) || is.null(experiments)) {
      experiments <- get_experiments(x)
    } else if (missing(range.lower) || missing(range.upper)) {
      warning("Please provide params 'range.lower' and 'range.upper'")
    }
    check_lengths(x, range.lower, range.upper)
    all.regions <- get_region_counts(x,
                                     regions,
                                     range.lower,
                                     range.upper,
                                     length = TRUE,
                                     transcript = TRUE,
                                     experiments = experiments)
  } else if (is.data.table(x)) {
    col.names <- c("experiment", "region.count", "region")
    mismatch <- !all(names(x) == col.names) || typeof(x[[1]]) != "character"||
      typeof(x[[2]]) != "double" || typeof(x[[3]]) != "character" || ncol(x) != 3
    if (mismatch) {
      stop("Please make sure that the data table is of the correct format.")
    } else if (!identical(unique(x$region), regions)) {
      stop("Please make sure that the data table only includes the 'UTR5','CDS', and 'UTR3' regions.")
    }
    all.regions <- x
  } else {
    stop("Please provide a ribo object.")
  }

  all.regions %>%
    group_by(.data$experiment) %>%
    summarize(sum = sum(.data$region.count)) -> total.counts

  total.counts %>%
    left_join(all.regions, by = "experiment") -> all.regions

  all.regions %>%
    mutate(fraction = .data$region.count/sum) %>%
    mutate(region = factor(.data$region, levels = c("UTR3", "CDS", "UTR5"))) %>%
    arrange(desc(.data$region)) %>%
    ggplot(aes(x = .data$experiment, y = .data$fraction, fill = .data$region)) +
    geom_col() +
    coord_flip() +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_fill_discrete(breaks = c("UTR5", "CDS", "UTR3")) +
    labs(title = title, x = "Experiment", y = "Fraction", fill = "Region")
}
