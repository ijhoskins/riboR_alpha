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
#' region.counts <- get_region_counts(ribo.object = sample,
#'                                    regions = regions,
#'                                    range.lower = 2,
#'                                    range.upper = 5,
#'                                    length = FALSE,
#'                                    transcript = TRUE,
#'                                    tidy = FALSE,
#'                                    experiments = experiments)
#'
#' @param ribo.object S3 "ribo" class object
#' @param range.lower Lower bound of the read length
#' @param range.upper Upper bound of the read length
#' @param length Option to condense the read lengths together, preserve the transcripts
#' @param transcript Option to condense the transcripts together, preserve the read lengths
#' @param tidy logical value denoting whether or not the user wants a tidy format
#' @param regions Specific region of interest
#' @param experiments List of experiment names
#' @return A data table of the region counts
#' @importFrom rhdf5 h5read
#' @importFrom data.table data.table
#' @importFrom tidyr gather
#' @export
get_region_counts <- function(ribo.object,
                              range.lower,
                              range.upper,
                              length = TRUE,
                              transcript = TRUE,
                              tidy = TRUE,
                              regions = c("UTR5", "UTR5J", "CDS", "UTR3J", "UTR3"),
                              experiments = get_experiments(ribo.object)) {

  range.info <- c(range.lower = range.lower,
                  range.upper = range.upper)

  regions <- check_rc_input(ribo.object,
                            regions,
                            range.info,
                            experiments)

  conditions <- c(transcript = transcript,
                  length     = length)

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

  values <- c("UTR5" = 1, "UTR5J" = 2, "CDS" = 3, "UTR3J" = 4, "UTR3" = 5)
  cols   <- unname(values[regions]) 
  total.regions <- length(regions)
  
  region.data <- determine_matrix_size(conditions,
                                       ref.length,
                                       total.regions,
                                       total.experiments,
                                       read.length.range)

  colnames(region.data) <- regions

  region.info <- list(ribo.min            = ribo.min,
                      range.info          = range.info,
                      ref.names           = ref.names,
                      total.experiments   = total.experiments,
                      matched.experiments = matched.experiments,
                      ref.length          = ref.length,
                      column.index        = cols,
                      conditions          = conditions)

  result <- fill_region(handle, region.info, region.data)
  
  if (tidy) {
    result <- setDT(gather(result, "region", "count", regions))
  }
  
  return(result)
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
    path <- paste("/experiments/", 
                  experiment,
                  "/region_counts/region_counts", sep ="")
    for (current.length in (range.lower:range.upper)) {
      offset <- ref.length * (current.length - ribo.min) + 1
      row.index    <- c(offset: (offset + ref.length - 1))
      region.count <- t(h5read(handle, 
                               path,
                               index = list(column.index, row.index)))
      
      fill.params <- list(index          = i,
                          conditions     = region.info[["conditions"]],
                          ref.length     = ref.length,
                          result         = region.data,
                          current.length = current.length,
                          range.info     = region.info[["range.info"]],
                          data           = region.count)

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
  regions <- check_regions(ribo.object, regions)

  range.lower <- range.info[["range.lower"]]
  range.upper <- range.info[["range.upper"]]
  
  check_lengths(ribo.object, range.lower, range.upper)
  check_experiments(ribo.object, experiments)
  return(regions)
}

check_regions <- function(ribo.object,
                          regions) {
  regions <- toupper(regions)
  region.options <- c("UTR5", "UTR5J", "CDS", "UTR3J", "UTR3")
  
  if (typeof(regions) != 'list' & typeof(regions) != "character") {
    stop("Please specify the regions as a single string, a vector, or a list.",
         call. = FALSE)
  }
  
  for (region in regions) {
    if (!region %in% region.options) {
      stop("'", region, "'", " is not an option. Please indicate the region with
           one of the following: 'UTR5', 'UTR5J, 'CDS', 'UTR3J', 'UTR3'",
           call. = FALSE)
    }
  }
 
  regions <- factor(regions, levels = c("UTR5", "UTR5J", "CDS", "UTR3J", "UTR3"), ordered = TRUE)
  regions <- as.vector(sort(regions))
  return(regions)
}



#' Plots the length distribution
#'
#' The function \code{\link{plot_length_distribution}} can take either a data.table
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
#' plot_length_distribution(x = sample,
#'                          region = "CDS",
#'                          range.lower = 2,
#'                          range.upper = 5,
#'                          experiments = experiments,
#'                          percentage = TRUE)
#'
#'
#' #data.table use case
#' #obtains the region counts at each individual read length, summed across every transcript
#' region.counts <- get_region_counts(ribo.object = sample,
#'                                    range.lower = 2,
#'                                    range.upper = 5,
#'                                    length = FALSE,
#'                                    tidy = TRUE,
#'                                    regions = regions,
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
#' @importFrom stats aggregate
#' @importFrom rlang .data sym
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
                           range.lower = range.lower,
                           range.upper = range.upper,
                           regions = region,
                           length = FALSE,
                           transcript = TRUE,
                           experiments = experiments)
  } else if (is.data.table(x)){
    col.names <- c("experiment", "length", "region", "count")
    mismatch <- !all(names(x) == col.names) || typeof(x[[1]]) != "character"||
    typeof(x[[2]]) != "integer"  || typeof(x[[3]]) != "character" ||
    typeof(x[[4]]) != "double" || ncol(x) != 4
    if (mismatch) {
      stop("Please make sure that the data table is of the correct format.")
    }
  } else {
    stop("Please make sure that param 'x' is either a data.table or a ribo object.")
  }

  y.axis <- "Count"
  y.value <- sym("count")

  if (percentage) {
    total <- aggregate(x$count, by = list(experiment= x$experiment), FUN=sum)
    x <- x %>%  left_join(total, by= "experiment") %>%
      mutate(fraction = .data$count/x) %>%
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
#'                                    regions = regions,
#'                                    range.lower = 2,
#'                                    range.upper = 5,
#'                                    tidy = TRUE,
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
#' @importFrom ggplot2 labs scale_fill_discrete element_blank geom_text position_stack
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
    all.regions <- get_region_counts(ribo.object = x,
                                     regions     = regions,
                                     range.lower = range.lower,
                                     range.upper = range.upper,
                                     length      = TRUE,
                                     transcript  = TRUE,
                                     experiments = experiments)
  } else if (is.data.table(x)) {
    col.names <- c("experiment", "region", "count")
    mismatch <- !all(names(x) == col.names) || typeof(x[[1]]) != "character"||
      typeof(x[[2]]) != "character" || typeof(x[[3]]) != "double" || ncol(x) != 3
    if (mismatch) {
      stop("Please make sure that the data table is of the correct format.")
    } else if (!identical(unique(x$region), regions)) {
      stop("Please make sure that the data table only includes the 'UTR5','CDS', and 'UTR3' regions.")
    }
    all.regions <- x
  } else {
    stop("Please provide a ribo object or a data.table of the correct format.")
  }
  
  #prepare data for visualization
  all.regions %>%
    group_by(.data$experiment) %>%
    summarize(sum = sum(.data$count)) -> total.counts

  total.counts %>%
    left_join(all.regions, by = "experiment") -> all.regions

  all.regions %>%
    mutate(percentage = round(100 * .data$count/sum, 1)) %>%
    mutate(region = factor(.data$region, levels = c("UTR3", "CDS", "UTR5"))) ->
    all.regions
  
  #text label of percentages only in the "CDS" region
  percentages <- replace(all.regions$percentage, all.regions$region != "CDS", "") 
  
  all.regions %>% 
    ggplot(aes(x    = .data$experiment, 
               y    = .data$percentage, 
               fill = .data$region)) +
      geom_col() +
      coord_flip() +
      geom_text(aes(x     = .data$experiment, 
                    y     = 50, 
                    label = percentages),
                size = 3) +
      theme_bw() +
      theme(plot.title   = element_text(hjust = 0.5),
            panel.border = element_blank(),
            panel.grid   = element_blank()) +
      scale_fill_discrete(breaks = c("UTR5", "CDS", "UTR3")) +
      labs(title = title,
           fill  = "Region",
           x     = "Experiment", 
           y     = "Percentage")
}