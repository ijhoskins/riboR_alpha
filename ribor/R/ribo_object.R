#' Creates an S3 object of class "ribo"
#'
#' \code{\link{ribo}} creates an S3 object. It creates a handle, extracts the root folder attributes,
#' and provides information about the reference transcript names and lengths
#'
#' This object is required as an argument for almost all of the functions in this package, and all of the
#' functions in this package can accept the returned object of this function. This object is not meant to 
#' be modified or changed by the user. It is meant to serve as an intermediary between the .ribo file and 
#' an R environment by creating an object that holds pertinent information.
#' 
#' The information stored in this object include the .ribo file handle, the list of experiments,
#' the format version, the reference model, the maximum read length, the minimum read length, the left span,
#' the right span, and other information about the transcript information.
#' 
#' Some of the subsequent function calls avoid direct usage of the information stored in the 
#' .ribo object to prevent any accidental error. However, certain variables in the returned object, such 
#' as the handle, are required to make use of the additional functionality in this package.
#' 
#' @param name The path to the .ribo file
#' @return Returns a list containing a handle to the HDF5 file,
#'        various attributes in the root folder, and
#'        information about the transcripts such as
#'        names and lengths
#' @importFrom rhdf5 H5Fopen h5readAttributes h5ls h5read
#' @importFrom hash hash
#' @export
ribo <- function(name){
  ribo.handle   <- H5Fopen(name)
  attributes <- h5readAttributes(ribo.handle, name = "/")

  transcript.names   <- h5read(ribo.handle&'reference',
                              name = "reference_names")
  
  transcript.lengths <- h5read(ribo.handle&'reference',
                              name = "reference_lengths")
  num.transcripts <- length(transcript.names)
  
  hash.value <- rep(list(c("offset" = 0, "length" = 0)), length = num.transcripts)
  names(hash.value) <- transcript.names
  hash.value[[1]]["length"] <-  transcript.lengths[[1]]
  
  for (i in 2:num.transcripts) {
    hash.value[[i]][["offset"]] <- hash.value[[i - 1]][["length"]] + hash.value[[i - 1]][["offset"]]
    hash.value[[i]][["length"]] <- transcript.lengths[[i]]
  }
  
  transcript.info <- hash(hash.value)
  length.offset <- hash.value[[num.transcripts]][["offset"]] +
                   hash.value[[num.transcripts]][["length"]]
  
  names(length.offset) <- NULL
  
  transcript.info <- hash(hash.value)
  
  ribo.contents <- list(handle             = ribo.handle,
                        experiments        = h5ls(ribo.handle&'experiments', recursive = FALSE)$name,
                        format.version     = attributes$format_version,
                        reference          = attributes$reference,
                        length.max         = attributes$length_max,
                        length.min         = attributes$length_min,
                        left.span          = attributes$left_span,
                        right.span         = attributes$right_span,
                        length.offset      = length.offset,
                        experiment.info    = get_content_info(ribo.handle),
                        transcript.info    = transcript.info)
  attr(ribo.contents, "class") <- "ribo"
  return(ribo.contents)
}

#' Printing "ribo" objects 
#' 
#' Print a ribo object
#' 
#' This method retrieves any user-relevant information about the ribo object 
#' and neatly prints it out. The information is extracted from the "ribo" class
#' object, turned into a data.table. The entries are then converted to a character,
#' customly formatted, and printed.
#' 
#' @param x object of class "ribo"
#' @param ... additional arguments that have no function yet 
#' @importFrom data.table as.data.table data.table
#' @importFrom hash keys
#' @export
print.ribo <- function(x, ...) {
  check_ribo(x)
  attributes <- c("format version", "reference", "min read length",
                  "max read length", "left span", "right span", "transcript count")
  left.span  <- x$left.span 
  right.span <- x$right.span 
  min.length <- x$length.min
  max.length <- x$length.max
  format.version <- x$format.version
  reference      <- x$reference 
  transcripts <- length(x$transcript.info)
  experiment.info <- x$experiment.info
  
  file.values <- c(format.version, reference, min.length, max.length,
                 left.span, right.span, transcripts)
  file.info   <- data.table("info" = attributes,
                            " " = file.values)
  
  name.width <- max(max(sapply(file.info, nchar)), max(sapply(names(file.info), nchar)))
  names(file.info) <- format(names(file.info), width = name.width, justify = "centre")
  file.info <- format(file.info, width = name.width, justify = "centre")
  
  #formatting the experiment.info
  experiment.info <- lapply(experiment.info, as.character)
  experiment.info <- as.data.table(experiment.info)
  name.size <- max(sapply(experiment.info[, -1], nchar))
  name.size <- max(name.size, max(sapply(names(experiment.info)[-1], nchar)))
  names(experiment.info)[-1] <- format(names(experiment.info)[-1], width = name.size, justify = "centre")
  experiment.info <- format(experiment.info, width = name.size, justify = "centre")
  experiment.info <- as.data.table(experiment.info)
  
  name.size <- max(sapply(experiment.info[, 1], nchar), nchar(names(experiment.info[1])))
  names(experiment.info)[1] <- format(names(experiment.info)[1], width = name.size, justify = "centre")
  experiment.info[, 1] <- format(experiment.info[, 1], width = name.size, justify = "centre")
               
  cat("General File Information:\n")
  print(data.table(file.info), row.names = FALSE, quote = FALSE)
  cat("\n")
  cat("Dataset Information:\n")
  print(data.table(experiment.info), row.names = FALSE, quote = FALSE)
  invisible(x)
}