#' Retrieves the coverage data for a given transcript
#'
#' The function \code{\link{get_coverage}} generates a data.table of coverage
#' data over the length of a given transcript.
#'
#' The function \code{\link{get_coverage}} first checks the experiments in the
#' 'experiment.list' parameter to see if they are present in the .ribo file.
#' It will then check these experiments for coverage data which is an optional
#' dataset. The function checks the coverage of one transcript at a time at
#' each read length from 'range.lower' to 'range.upper', inclusive. However,
#' the parameter 'length' allows the user to obtain the coverage
#' information of a transcript across the range of read lengths indicated by
#' 'range.lower' and 'range.upper'.
#' @examples
#' #generate the ribo object
#' file.path <- system.file("extdata", "sample.ribo", package = "ribor")
#' sample <- ribo(file.path)
#'
#' #get the experiments of interest that also contain coverage data
#' experiments <- c("Hela_1", "Hela_2", "Hela_3", "WT_1")
#'
#' #the ribo file contains a transcript named 'MYC'
#' coverage.data <- get_coverage(ribo.object = sample,
#'                               name = "MYC",
#'                               range.lower = 2,
#'                               range.upper = 5,
#'                               length = TRUE,
#'                               experiments = experiments)
#'
#' @param ribo.object S3 "ribo" class object
#' @param name Name of the transcript
#' @param range.lower Lower bound of the read length
#' @param range.upper Upper bound of the read length
#' @param length Logical value that denotes if the coverage should be summed across read lengths
#' @param experiments List of experiments to obtain coverage information on
#' @return A data table containing the coverage data
#' @seealso \code{\link{ribo}} to generate the necessary ribo.object parameter
#' @importFrom rhdf5 h5read
#' @export
get_coverage <- function(ribo.object,
                         name,
                         range.lower,
                         range.upper,
                         length = TRUE,
                         experiments = get_experiments(ribo.object)) {
  #perform checks on the parameters
  check_lengths(ribo.object, range.lower, range.upper)
  check_experiments(ribo.object, experiments)
  ribo.experiments    <- get_experiments(ribo.object)

  #generate list of experiments also present in the ribo file that have coverage data
  filter.experiments  <- intersect(experiments, ribo.experiments)
  coverage.list       <- check_coverage(ribo.object, filter.experiments)
  matched.experiments <- intersect(filter.experiments, coverage.list)
  total.experiments   <- length(matched.experiments)

  #generate offsets
  current.offset      <- ribo.object$transcript.offset[[name]]
  length.offset       <- ribo.object$length.offset
  transcript.length   <- ribo.object$transcript.lengths[name]
  min.length          <- get_read_lengths(ribo.object)[1]
  read.range          <- range.upper - range.lower + 1

  #create matrix of correct size based on param transcript
  result <- matrix()
  if (length) {
    result <- matrix(0L,
                     nrow = 1 * total.experiments,
                     ncol = transcript.length)
  } else {
    result <- matrix(nrow = read.range * total.experiments,
                     ncol = transcript.length)
  }
  colnames(result) <- c(1:transcript.length)

  #fills in the entries of the matrix
  for (i in 1:total.experiments) {
    #for the given experiment, get its coverage
    experiment <- matched.experiments[i]
    path <- paste("/experiments/", experiment, "/coverage/coverage", sep ="")

    #get the coverage information for the specific transcript at each read length
    for (current.length in range.lower:range.upper) {
      correct.length <- 1 + (current.length - min.length) * length.offset
      coverage.start <- correct.length + current.offset
      coverage.stop  <- coverage.start + transcript.length - 1
      coverage<- t(h5read(ribo.object$handle,
                          path,
                          index = list(coverage.start:coverage.stop)))
      coverage <- as.integer(coverage)
      # compute the correct offset to store the coverage information in each case
      if (length) {
        result[i, ] <- result[i, ] + coverage
      } else {
        current.experiment      <- (i - 1) * read.range
        current.read            <- current.length - range.lower
        current.index           <- current.experiment + current.read + 1
        result[current.index, ] <- coverage
      }
    }
  }
  return(create_datatable(length,
                          range.lower,
                          range.upper,
                          matched.experiments,
                          result))
}


create_datatable <- function (length,
                              range.lower,
                              range.upper,
                              matched.experiments,
                              matrix) {
  # Given a matrix of coverage data, create_datatable generates the correct
  # data.table based on a set of parameters
  #
  # Returns:
  # Data table that wraps the matrix with the correct and appropriate labels

  matched.size <- length(matched.experiments)
  if (length) {
    return (data.table(experiment = matched.experiments,
                       matrix))
  }
  range <- range.upper - range.lower + 1
  return (data.table(experiment  = rep(matched.experiments, each = range),
                     read.length = rep(c(range.lower:range.upper), matched.size),
                     matrix))
}


check_coverage <- function(ribo.object, experiments) {
  # helper function that both generates a list of experiments with coverage data using get_info
  # and produces a warning message for each experiment(provided by the user) that does not have
  # coverage data
  #
  # Args:
  # ribo.object: S3 object of class "ribo"
  # experiment.list: list of experiments inputted by the user
  #
  # Returns:
  # A list of experiments in the ribo.object that have coverage data

  #obtain the coverage data
  table <- get_info(ribo.object)[["contents"]][, c("names", "coverage")]
  has.coverage <- table[table$coverage == TRUE, ]
  has.coverage <- has.coverage$names

  #find the experiments in the experiment.list that do not have coverage and print warnings
  check <- setdiff(experiments, has.coverage)
  if (length(check)) {
    for (experiment in check) {
      warning("'", experiment, "'", " did not have coverage data.")
    }
    warning("Param experiment.list contains experiments that did not have coverage data.
            The return value ignores these experiments.")
  }

  #return a list of experiments with coverage
  return(has.coverage)
}
