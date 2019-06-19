check_ribo <- function(ribo.object, stop = TRUE) {
  # Helper method that checks the internal contents and class of a parameter
  # check_ribo takes in an object and checks for proper contents and class
  # Args:
  # ribo.object An S3 object of class "ribo"
  #
  # Return:
  # None

  expected_names <- c("handle", "experiments", "format.version", "reference",
                      "length.max", "length.min", "left.span", "right.span",
                      "length.offset", "transcript.offset", "transcript.names",
                      "transcript.lengths")
  if(!identical(names(ribo.object), expected_names) || class(ribo.object) != "ribo") {
    if (stop) {
      stop("Param ribo.object should be of class ribo.", call. = FALSE)
    }
    return(FALSE)
  }
  return(TRUE)
}

check_lengths <- function(ribo.object, range.lower, range.upper) {
  # Helper method that checks for correct lengths
  #
  # check_lengths directly reads the .ribo file for its lowest and highest read
  # length and compares it to the corresponding parameters
  #
  # Args:
  # ribo.object S3 object of class "ribo"
  # range.lower lowest read length
  # range.upper highest read length
  #
  # Return:
  # none
  min.length <- get_attributes(ribo.object)$length_min
  max.length <- get_attributes(ribo.object)$length_max
  if ((range.lower < min.length | range.lower > range.upper)){
    stop("Param range.lower must be greater than or equal to the minimum
          length and less than range.upper.", call. = FALSE)
  } else if ((range.upper > max.length | range.upper < range.lower)) {
    stop("Param range.upper must be less than or equal to the maximum length
          and greater than or equal to range.lower.", call. = FALSE)
  }
}


check_experiments <- function(ribo.object, experiments) {
  # Helper method that checks if the user-given experiments
  # are present in the current ribo file
  # Args:
  # ribo.object - S3 object of class ribo, contains the handle to the file
  # Return:
  # None
  ribo.experiments <- get_experiments(ribo.object)
  matched.experiments <- intersect(experiments, ribo.experiments)

  #deals with missing experiments
  missing <- FALSE
  check <- setdiff(matched.experiments, experiments)
  if (length(check)) {
    for (experiment in check) {
      missing <- TRUE
      warning("'", experiment, "'", " was not found.", call. = FALSE)
    }
    warning("Param 'experiments' contained experiments that were not found.
            The returned data table ignores these experiments.", call. = FALSE)
  }
}
