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
#' @export
ribo <- function(name){
  ribo.handle   <- H5Fopen(name)
  attributes <- h5readAttributes(ribo.handle, name = "/")

  transcript.names   <- h5read(ribo.handle&'reference',
                              name = "reference_names")
  transcript.lengths <- h5read(ribo.handle&'reference',
                              name = "reference_lengths")
  names(transcript.lengths) <- transcript.names
  total.ref <- length(transcript.names)
  transcript.offset  <- vector("list", total.ref)
  transcript.offset[[1]] <- 0
  names(transcript.offset)[1] <- transcript.names[1]
  for (i in 2:total.ref) {
    transcript.offset[[i]] <- transcript.offset[[i-1]] + transcript.lengths[i-1]
    names(transcript.offset)[i] <- transcript.names[i]
  }
  length.offset <- transcript.offset[[total.ref]] +
                   transcript.lengths[total.ref]

  ribo.contents <- list(handle             = ribo.handle,
                        experiments        = h5ls(ribo.handle&'experiments', recursive = FALSE)$name,
                        format.version     = attributes$format_version,
                        reference          = attributes$reference,
                        length.max         = attributes$length_max,
                        length.min         = attributes$length_min,
                        left.span          = attributes$left_span,
                        right.span         = attributes$right_span,
                        length.offset      = length.offset,
                        transcript.offset  = transcript.offset,
                        transcript.names   = transcript.names,
                        transcript.lengths = transcript.lengths)
  attr(ribo.contents, "class") <- "ribo"
  return(ribo.contents)
}