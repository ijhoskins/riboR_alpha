#' Retrieves the metagene data from a .ribo file
#'
#' The function \code{\link{get_metagene}} returns a data table that provides
#' the coverage at the positions surrounding the metagene start or stop site.
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
#' The param 'transcripts' condenses the transcripts together.
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
#' @param ribo.object S3 "ribo" class object
#' @param site "start" or "stop" site coverage
#' @param range.lower Lower bound of the read length, inclusive
#' @param range.upper Upper bound of the read length, inclusive
#' @param transcript Option to condense the transcripts together, preserving information at each read length
#' @param length Option to condense the read lengths together, preserving information at each transcript
#' @param experiments List of experiment names
#' @return A data.table of the metagene information
#' @examples
#'
#' #generate the ribo object by providing the file.path to the ribo file
#' file.path <- system.file("extdata", "sample.ribo", package = "ribor")
#' sample <- ribo(file.path)
#'
#'
#' #extract the total metagene information for all experiments
#' #across the read lengths and transcripts of the start site
#' #from read length 2 to 5
#' metagene_info <- get_metagene(ribo.object = sample,
#'                               site = "start",
#'                               range.lower = 2,
#'                               range.upper = 5,
#'                               length = TRUE,
#'                               transcript = TRUE,
#'                               experiments = get_experiments(sample))
#'
#'
#' #Note that length, transcript, and experiments in this case are the
#' #default values and can be left out. The following generates the same output.
#'
#' \donttest{
#' metagene_info <- get_metagene(ribo.object = sample,
#'                               site = "start",
#'                               range.lower = 2,
#'                               range.upper = 5)
#' }
#'
#'
#' @seealso
#' \code{\link{ribo}} to generate the necessary 'ribo' class object,
#' \code{\link{plot_metagene}} to visualize the metagene data,
#' \code{\link{get_tidy_metagene}} to obtain tidy metagene data under certain conditions
#' @importFrom rhdf5 h5read
#' @importFrom data.table data.table
#' @export
get_metagene <- function(ribo.object,
                         site,
                         range.lower,
                         range.upper,
                         length = TRUE,
                         transcript = TRUE,
                         experiments = get_experiments(ribo.object)) {

  range.info <- c(range.lower = range.lower,
                  range.upper = range.upper)

  site <- tolower(site)

  check_metagene_input(ribo.object,
                       site,
                       range.info,
                       experiments)

  conditions <-  c(length = length,
                   transcript = transcript)

  #gather information to use in filling and labeling final data table
  ribo.experiments <- get_experiments(ribo.object)
  handle <- ribo.object$handle
  matched.experiments <- intersect(experiments, ribo.experiments)
  total.experiments   <- length(matched.experiments)

  ref.names    <- get_reference_names(ribo.object)
  ref.length   <- length(get_reference_names(ribo.object))
  read.lengths <- get_read_lengths(ribo.object)
  ribo.min <- read.lengths[1]
  read.length.range   <- range.upper - range.lower + 1


  #compute the number of columns
  metagene.radius     <- h5readAttributes(handle, "/")[["metagene_radius"]]
  total.indices       <- 2 * metagene.radius + 1
  columns             <- c(1: total.indices)


  #get matrix of appropriate size
  result <- determine_matrix_size(conditions,
                                  ref.length,
                                  total.indices,
                                  total.experiments,
                                  read.length.range)

  colnames(result) <- c(-metagene.radius:metagene.radius)

  #loop over the experiments and fill in the matrix
  for (i in 1:total.experiments) {
      experiment <- matched.experiments[i]
      dataset.name <- paste(site, "_site_coverage", sep = "")
      datagroup.path <- paste("/experiments/", experiment, "/metagene/", sep ="")
      path <- paste(datagroup.path, dataset.name, sep = "")

      #read in all of the columns of a transcript for a given read length
      for (current.length in (range.lower:range.upper)) {
        offset <- ref.length * (current.length - ribo.min) + 1
        rows <- c(offset : (offset + ref.length - 1))
        metagene.data <- t(h5read(handle,
                                  path,
                                  index = list(columns,rows)))
        fill.params <-  fill.params <- list(index = i,
                                            conditions = conditions,
                                            ref.length = ref.length,
                                            result = result,
                                            current.length = current.length,
                                            range.info = range.info,
                                            data = metagene.data)

        result <- fill_matrix(fill.params)
      }
  }
  return(make_datatable(ref.names,
                        transcript,
                        length,
                        matched.experiments,
                        range.lower,
                        range.upper,
                        result))
}

#' Retrieves the metagene data in a tidy format
#'
#' The function \code{\link{get_tidy_metagene}} provides the user with a tidy data format for easier
#' data cleaning and manipulation. In providing this functionality, the user must length the
#' transcripts together and is only provided the option to length the read lengths together.
#'
#' The dimensions of the returned data table depend on the parameters
#' range.lower, range.upper, and length.
#'
#' The param 'length' condenses the read lengths together.
#' When length is true, then the resulting data table prints out one row
#' for each experiment. This provides a tidy format of the metagene information
#' across all transcripts and all read lengths in a given experiment. Each row
#' in the data table represents the total metagene coverage count of a given experiment
#' at a given position.
#'
#' When the param  'length' is false, then the resulting data table prints out the
#' metagene coverage count at each position of the metagene radius for each read length.
#' This provides a tidy format of the metagene information across the transcripts, preserving
#' the metagene coverage count at each read length.
#'
#' @param ribo.object S3 "ribo" class object
#' @param site "start" or "stop" site coverage
#' @param range.lower Lower bound of the read length
#' @param range.upper Upper bound of the read length
#' @param length Option to condense the read lengths together
#' @param experiments List of experiment names
#' @return
#' A tidy data table of the metagene information
#' @examples
#' #generate the ribo object by loading in a ribo function and calling the \code{\link{ribo}} function
#' file.path <- system.file("extdata", "sample.ribo", package = "ribor")
#' sample <- ribo(file.path)
#'
#' #extract the total metagene information in a tidy format
#' #for all experiments across the read lengths and transcripts
#' #of the start site from read length 2 to 5
#'
#' metagene_info <- get_tidy_metagene(ribo.object = sample,
#'                                    site = "start",
#'                                    range.lower = 2,
#'                                    range.upper = 5,
#'                                    length = TRUE,
#'                                    experiments = get_experiments(sample))
#'
#' #Note that length and experiments in this case are the
#' #default values and can be left out. The following generates the same output.
#' \donttest{
#' metagene_info <- get_tidy_metagene(ribo.object = sample,
#'                                    site = "start",
#'                                    range.lower = 2,
#'                                    range.upper = 5)
#' }
#'
#' @seealso
#' \code{\link{ribo}} to generate the necessary 'ribo' class object.
#' \code{\link{plot_metagene}} to visualize the metagene data,
#' \code{\link{get_metagene}} to obtain tidy metagene data under certain conditions
#' @importFrom rhdf5 h5read
#' @importFrom data.table data.table setDT
#' @export
get_tidy_metagene <- function(ribo.object,
                              site,
                              range.lower,
                              range.upper,
                              length = TRUE,
                              experiments = get_experiments(ribo.object)) {
  site <- tolower(site)
  result <- get_metagene(ribo.object,
                         site,
                         range.lower,
                         range.upper,
                         length,
                         transcript = TRUE,
                         experiments = experiments)
  metagene.radius <- as.integer((ncol(result) - 2)/2)
  tidy.data <- gather(result, key = "position", value = "count", c(as.character(-metagene.radius:metagene.radius)))
  tidy.data$position <- as.integer(tidy.data$position)
  tidy.data$count    <- as.integer(tidy.data$count)
  tidy.data %>% 
    left_join(ribo.object$experiment.info[, c("experiment", "total.reads")],
              by = "experiment") -> tidy.result
  return(setDT(tidy.result))
}


check_metagene_input <- function(ribo.object,
                                 site,
                                 range.info,
                                 experiments) {
  #check_metagene_input is a helper function that checks the validity of
  #the metagene function parameters

  #check param validity
  if (site != "start" & site != "stop") {
    stop("Please type 'start' or 'stop' to indicate the 'site' parameter value.")
  }

  range.lower <- range.info[["range.lower"]]
  range.upper <- range.info[["range.upper"]]

  check_ribo(ribo.object)
  check_lengths(ribo.object, range.lower, range.upper)
  check_experiments(ribo.object, experiments)
}

#' Plots the metagene coverage data
#'
#' The function \code{\link{plot_metagene}} plots the metagene site coverage,
#' separating by experiment.
#'
#' If a data.table is provided as param 'x', then the only additional parameter
#' is the optional title' parameter for the generated plot. If a ribo.object is
#' provided as param 'x', the rest of the parameters listed are necessary.
#'
#' When given a ribo class object, the \code{\link{plot_metagene}} function
#' generates a data.table by calling the \code{\link{get_tidy_metagene}}
#' function, so the run times in this case will be mostly comprised of a call
#' to the \code{\link{get_metagene}} function.
#'
#' This function uses ggplot in its underlying implementation.
#'
#' @param x ribo.object or data.table
#' @param site "start" or "stop" site
#' @param range.lower lower bound of the read length, inclusive
#' @param range.upper upper bound of the read length, inclusive
#' @param experiments list of experiments
#' @param normalize When TRUE, normalizes the data by the total reads.
#' @param title title of the generated plot
#' @param tick x-axis labeling increment
#' @examples
#' #a potential use case is to directly pass in the ribo object file as param 'x'
#'
#' #generate the ribo object to directly use
#' file.path <- system.file("extdata", "sample.ribo", package = "ribor")
#' sample <- ribo(file.path)
#'
#' #specify experiments of interest
#' experiments <- c("HeLa_1", "HeLa_2", "WT_1")
#'
#' #plot the metagene start site coverage for all experiments in 'sample.ribo'
#' #from read length 2 to 5
#' plot_metagene(x = sample,
#'               site = "start",
#'               range.lower = 2,
#'               range.upper = 5,
#'               experiments = experiments)
#'
#' #Note that the site, range.lower, range.upper, and experiments are only
#' #necessary if a ribo object is being passed in as param 'x'. If a ribo
#' #object is passed in, then the param 'experiments' will be set to all of
#' #the experiments by default.
#'
#' #If a data.table is passed in, then the plot_metagene function
#' #does not need any other information. All of the elements of the data.table
#' #will be used, assuming that it contains the same column names and number of
#' #columns as the output from get_tidy_metagene()
#'
#' #gets the metagene start site coverage from read length 2 to 5
#' #note that the data must be summed across transcripts and read lengths
#' #for the plot_metagene function
#' data <- get_tidy_metagene(sample,
#'                           site = "start",
#'                           range.lower = 2,
#'                           range.upper = 5)
#'
#' #plot the metagene data
#' plot_metagene(data)
#'
#'
#'
#' @importFrom data.table is.data.table
#' @importFrom dplyr left_join mutate
#' @importFrom ggplot2 ggplot geom_line theme_bw ggtitle aes
#' @importFrom ggplot2 element_text theme labs scale_x_continuous
#' @importFrom rlang .data sym
#' @importFrom tidyr gather
#' @export
plot_metagene <- function(x,
                          site,
                          range.lower,
                          range.upper,
                          experiments,
                          normalize = FALSE,
                          title = "Metagene Site Coverage",
                          tick = 10) {
  is.ribo <- check_ribo(x, stop = FALSE)
  
  #x is a ribo object
  if (is.ribo) {
    missing.ranges <- missing(range.lower) || missing(range.upper)
    if (missing.ranges) {
      stop("Please indicate the 'range.lower' and 'range.upper' parameters.")
    } else if (missing(site)) {
      stop("Please indicate the 'site' parameter with either 'start' or 'stop'")
    } else if (missing(experiments)) {
      experiments = get_experiments(x)
    }
    x <- get_tidy_metagene(x,
                           site,
                           range.lower,
                           range.upper,
                           length = TRUE,
                           experiments = experiments)
  } else if (is.data.table(x)){
    metagene.radius <- (ncol(x) - 2)/2
    col.names       <- c("experiment", "position", "count", "total.reads")

    mismatch <- !all(names(x) == col.names) || (typeof(x[[1]]) != "character")
    col <- 2
    while (!mismatch && col < ncol(x)) {
      mismatch <- typeof(x[[col]]) != "integer"
      col <- col + 1
    }
    if (mismatch) {
      stop("Data table is not of the correct form.")
    }
  } else { #not a data table
    stop("Please make sure that param 'x' is either a data.table or a ribo object.")
  }
  
  y.value <- sym("count")
  y.label <- "Count"
  if (normalize) {
    x %>%
      mutate(normalize = (.data$count/.data$total.reads) * 1000) -> x
    y.value <- sym("normalize")
    y.label <- "Normalized Count"
  }
  
  metagene.radius <- max(x$position)

  ggplot(x,
         aes(x     = .data$position,
             y     = !!y.value,
             color = .data$experiment)) +
    scale_x_continuous(breaks = seq(-metagene.radius, metagene.radius, tick)) +
    geom_line() +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5)) +
    labs(title = title, x = "Position", y = y.label, color = "Experiment")
}