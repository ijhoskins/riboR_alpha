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

  for (i in 2:total.ref) {
    transcript.offset[[i]] <- transcript.offset[[i-1]] + transcript.lengths[[i-1]]
  }
  names(transcript.offset) <- transcript.names
  length.offset <- transcript.offset[[total.ref]] +
                   transcript.lengths[total.ref]
  #creates the separate lists for reads, coverage, rna.seq, and metadata
  #to eventually put in a data frame
  reads.list    <- list()
  coverage.list <- list()
  rna.seq.list  <- list()
  metadata.list <- list()
  
  #ls function provides information about the contents of each experiment
  ls <- h5ls(ribo.handle)
  experiments        = h5ls(ribo.handle&'experiments', recursive = FALSE)$name
  #loop over all of the experiments
  for (experiment in experiments) {
    #gathers information on the number of reads for each experiment by looking at
    #the attributes
    name           <- paste("/experiments/", experiment, sep = "")
    attribute      <- h5readAttributes(ribo.handle, name)
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
  
  experiments.info       <- data.table(experiment    = experiments,
                                       reads    = reads.list,
                                       coverage = coverage.list,
                                       rna.seq  = rna.seq.list,
                                       metadata = metadata.list)
  
  ribo.contents <- list(handle             = ribo.handle,
                        experiments        = experiments,
                        format.version     = attributes$format_version,
                        reference          = attributes$reference,
                        length.max         = attributes$length_max,
                        length.min         = attributes$length_min,
                        left.span          = attributes$left_span,
                        right.span         = attributes$right_span,
                        length.offset      = length.offset,
                        experiment.info    = experiments.info,
                        transcript.offset  = transcript.offset,
                        transcript.names   = transcript.names,
                        transcript.lengths = transcript.lengths)
  attr(ribo.contents, "class") <- "ribo"
  return(ribo.contents)
}

#' @export
print.ribo <- function(ribo.object) {
  values <- "File Information"
  attributes <- c("format version", "reference", "min read length",
                  "max read length", "left span", "right span", "transcript count")
  left.span  <- ribo.object$left.span 
  right.span <- ribo.object$right.span 
  min.length <- ribo.object$length.min
  max.length <- ribo.object$length.max
  format.version <- ribo.object$format.version
  reference      <- ribo.object$reference 
  transcripts <- length(ribo.object$transcript.names)
  experiment.info <- ribo.object$experiment.info
  
  file.values <- c(format.version, reference, min.length, max.length,
                 left.span, right.span, transcripts)
  file.info   <- data.table("---info--- " = attributes,
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
  invisible(ribo.object)
}

