#' Get information about the .ribo file
#'
#' The function \code{\link{get_info}} provides information on the attributes, metadata,
#' and datasets of the ribo file.
#'
#' The \code{\link{get_info}} first provides information on the format version, left_span, right_span,
#' longest read length, shortest read length, metagene_radius, and reference model. The last element of the
#' returned list contains the information about the presence of coverage and RNA-seq data which are
#' optional datasets to include in a .ribo file.
#'
#' @param ribo.object ribo.object is an S3 object of class "ribo"
#' @examples
#' #generate the ribo object
#' file.path <- system.file("extdata", "sample.ribo", package = "ribor")
#' sample <- ribo(file.path)
#'
#' #retrieve information
#' get_info(sample)
#'
#' @return Returns a list of 8 elements providing the following information in order of
#' format version, left_span, length_max, length_min, metagene_radius, reference,
#' right_span,and a data table that provides information about each experiment, including
#' names, number of reads, metadata, and its internal datasets
#'
#' @seealso \code{\link{ribo}} to generate the necessary ribo.object parameter
#'
#' @importFrom rhdf5 h5ls h5readAttributes
#' @importFrom data.table data.table
#' @export
get_info <- function(ribo.object) {
  check_ribo(ribo.object)
  handle   <- ribo.object$handle

  #retrieve an experiment list
  exp.list <- get_experiments(ribo.object)
  result <- get_attributes(ribo.object)
  if ("time" %in% names(result)) {
    result <- result[ -which(names(result) == "time")]
  }

  #creates the separate lists for reads, coverage, rna.seq, and metadata
  #to eventually put in a data frame
  reads.list    <- list()
  coverage.list <- list()
  rna.seq.list  <- list()
  metadata.list <- list()

  #ls function provides information about the contents of each experiment
  ls <- h5ls(handle)

  #loop over all of the experiments
  for (experiment in exp.list) {
    #gathers information on the number of reads for each experiment by looking at
    #the attributes
    name           <- paste("/experiments/", experiment, sep = "")
    attribute      <- h5readAttributes(handle, name)
    reads.list     <- c(reads.list, attribute[["total_reads"]])


    #creates separate logical lists to denote the presence of
    #reads, coverage, RNA-seq, metadata
    has.metadata   <- ("metadata" %in% names(attribute))
    metadata.list  <- c(metadata.list, has.metadata)

    group.contents <- ls[ls$group == name,]
    group.names    <- group.contents$name

    has.coverage   <- ("coverage" %in% group.names)
    coverage.list  <- c(coverage.list, has.coverage)

    has.rna.seq    <- ("rnaseq" %in% group.names)
    rna.seq.list   <- c(rna.seq.list, has.rna.seq)
  }

  #contructs the data.table containing logical values to denote the contents
  experiments.info       <- data.table(names    = exp.list,
                                       reads    = reads.list,
                                       coverage = coverage.list,
                                       rna.seq  = rna.seq.list,
                                       metadata = metadata.list)

  experiments                   <- list(experiments.info)
  result                        <- c(result, contents = experiments)
  return(result)
}


#' Retrieves the metadata of an experiment
#'
#' \code{\link{get_metadata}} provides information on the cell_line, digestion_duration,
#' digestion_enzyme, and a website link to the source if the experiment is found in the .ribo
#' file. If the experiment is not found, then the attributes of the root .ribo file is returned instead.
#'
#' @param ribo.object S3 class of object ribo
#' @param name The name of the experiment
#' @examples
#' #ribo object use case
#' #generate the ribo object
#' file.path <- system.file("extdata", "sample.ribo", package = "ribor")
#' sample <- ribo(file.path)
#'
#' #the ribo file contains an experiment named 'Hela_1'
#' get_metadata(sample, "Hela_1")
#'
#' @return
#' If the experiment is valid, a list of 4 elements providing
#' information about the cell_line, digestion_duration, digestion_enzyme, and
#' source link is returned.
#'
#' If the name is not found, then a list of attributes describing the root file
#' is provided.
#' @seealso \code{\link{ribo}} to generate the necessary ribo.object parameter
#' @importFrom rhdf5 h5readAttributes
#' @importFrom yaml read_yaml
#' @export
get_metadata <- function(ribo.object, name) {
  #get_experiments also checks if ribo.object is a proper "ribo" object
  exp.list <- get_experiments(ribo.object)
  handle <- ribo.object$handle

  #check if the name is in the list
  if (!(name %in% exp.list)) {
    warning("Experiment was not found. Returning root file attributes.")
    attribute <- h5readAttributes(handle, "/")
    if ("time" %in% names(result)) {
      result <- result[ -which(names(result) == "time")]
    }
    return (attribute)
  }

  #create the experiment path and get its attributes
  path <- paste("experiments/", name, sep = "")
  attribute <- h5readAttributes(handle, path)

  #check for metadata
  if ("metadata" %in% names(attribute)) {
    read_yaml(text = attribute[["metadata"]])
  } else {
    warning("'", name, "'", " does not have metadata. Returning an empty list")
    return(list())
  }
}

#' Provides a list of experiments from a .ribo file
#'
#' The function \code{\link{get_experiments}} provides a list of experiment names in the .ribo file.
#'
#' \code{\link{get_experiments}} returns a list of strings denoting the experiments. It obtains this
#' by reading directly from the .ribo file through the handle of the 'ribo.object' parameter. To generate
#' the param 'ribo.object', call the \code{\link{ribo}} function and provide the path to the .ribo file of interest.
#'
#' The user can then choose to create a subset from this list for any specific experiments of interest
#' for later function calls. Many functions that have the param 'experiment.list'
#'  call \code{\link{get_experiments}} to generate a default list of all experiments in the
#' .ribo file.
#' @examples
#' #generate the ribo object
#' file.path <- system.file("extdata", "sample.ribo", package = "ribor")
#' sample <- ribo(file.path)
#'
#' #get a list of the experiments
#' get_experiments(sample)
#'
#' @seealso \code{\link{ribo}} to generate the necessary ribo.object parameter
#' @param ribo.object S3 object of class "ribo"
#' @return A list of the experiment names
#' @importFrom rhdf5 h5ls
#' @export
get_experiments <- function(ribo.object) {
  check_ribo(ribo.object)
  result <- h5ls(ribo.object$handle)
  result <- result[result$group == "/experiments", ]
  return(result$name)
}
