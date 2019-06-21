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
  
  experiment.info <- list(ribo.object$experiment.info)
  result <- c(result, experiment.info = experiment.info)
  return(result)
}

print_info <- function(ribo.object) {
  info <- get_info(ribo.object)
  info.table <- data.table("--name--" = names(info[-length(info)]),
                           " " = c(as.character(info[-length(info)])))
  name.width <- max(max(sapply(info.table, nchar)), max(sapply(names(info.table), nchar)))
  names(info.table) <- format(names(info.table), width = name.width, justify = "centre")
  info.table <- format(info.table, width = name.width, justify = "left")
  cat("General Info:\n")
  print(data.table(info.table), row.names = FALSE, quote = FALSE)
  cat("\n")
  cat("Experiment Contents and Information:\n")
  contents <- info[[length(info)]]
  names(contents) <- format
  print(info[[length(info)]], row.names = FALSE, justify = "right")
}



#' Retrieves the metadata of an experiment
#'
#' \code{\link{get_metadata}} provides information on all of the user-inputted
#' metadata of an experiment. If the experiment is not found, then the 
#' attributes of the root .ribo file is returned instead.
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
#' If the experiment is valid, a list of elements providing all of the metadata of the  
#' experiment.
#'
#' If the name is not found, then a list of attributes describing the root file
#' is provided.
#' @seealso \code{\link{ribo}} to generate the necessary ribo.object parameter
#' @importFrom rhdf5 h5readAttributes
#' @importFrom yaml read_yaml
#' @export
get_metadata <- function(ribo.object, name = NULL) {
  #get_experiments also checks if ribo.object is a proper "ribo" object
  exp.list <- get_experiments(ribo.object)
  handle <- ribo.object$handle

  if (is.null(name)) {
    attributes <- get_attributes(ribo.object)
    return(data.table(attribute = names(attributes), value = attributes))
  } else if (!(name %in% exp.list)) {
    warning("'", name, "'", " is not a valid experiment name.",
            " Returned value is the root ribo file attributes.",
            call. = FALSE)
    attributes <- get_attributes(ribo.object)
    return(data.table(attribute = names(attributes), value = attributes))
  }
  
  #create the experiment path and get its attributes
  path <- paste("experiments/", name, sep = "")
  attribute <- h5readAttributes(handle, path)
  
  #check for metadata
  if ("metadata" %in% names(attribute)) {
    raw.result <- read_yaml(text = attribute[["metadata"]])
    result <- data.table(info = names(raw.result), " " = raw.result)
    
    filter <- result[result$info != "link", ]
    link <- unlist(result[result$info == "link", ][[2]])
    print(filter, row.names = FALSE)
    cat("link: \n")
    cat(link)
    invisible(raw.result)
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
