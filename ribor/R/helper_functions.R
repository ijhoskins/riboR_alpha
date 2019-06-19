get_reference_names <- function(ribo.object){
  # Retrieves the reference transcript names 
  check_ribo(ribo.object)
  return(h5read(ribo.object$handle&'reference', name = "reference_names"))
}


get_reference_lengths <- function(ribo.object){
  # Retrieves the reference transcript lengths 
  check_ribo(ribo.object)
  row.names <- h5read(ribo.object$handle&'reference', name = "reference_names")
  lengths   <- h5read(ribo.object$handle&'reference', name = "reference_lengths")
  return(data.table(transcript = row.names, length = lengths))
}


get_attributes <- function(ribo.object) {
  # Retrieves the attributes of the ribo.object 
  handle  <- ribo.object$handle
  result  <- list()
  result  <- h5readAttributes(handle, "/")

  #time does not serve much use
  result <- result[ -which(names(result) == "time")]
  return(result)
}

get_read_lengths <- function(ribo.object) {
  # Retrieves the minimum and maximum read lengths
  #
  # get_read_lengths finds the minimum and maximum read lengths of the .ribo file
  attributes <- get_attributes(ribo.object)
  result <- c(attributes$length_min, attributes$length_max)
  return(result)
}

fill_matrix <- function(info) {
  # helper method of the get_region_counts function that fills in the matrix
  # Returns:
  # result - an updated version of the matrix that has been filled in
  length <- info[["conditions"]][["length"]]
  transcript <- info[["conditions"]][["transcript"]]
  ref.length <- info[["ref.length"]]
  data <- info[["data"]]
  index <- info[["index"]]
  result <- info[["result"]]
  range.lower <- info[["range.info"]][["range.lower"]]
  range.upper <- info[["range.info"]][["range.upper"]]
  current.length <- info[["current.length"]]

  #sums across each length
  if (length & !transcript) {
    row.start <- (index - 1) * ref.length + 1
    row.stop  <-  index * ref.length
    #add the current matrix to the result
    result[row.start:row.stop, ] <- result[row.start:row.stop, ] + data
  } else if (transcript) {
    temp <- colSums(data)
    if (length) {
      #add to current result
      result[index, ] <- result[index, ] + temp
    } else {
      #compute the row offset
      row.start <- (index - 1) * (range.upper - range.lower + 1)
      row.start <- row.start + 1
      row.start <- row.start + (current.length - range.lower)

      #fill in matrix
      result[row.start, ] <- temp
    }
  } else {
    #compute the row offset
    row.start <- (index - 1) * ref.length * (range.upper - range.lower + 1)
    add <- (current.length - range.lower) * ref.length
    row.start <- row.start + 1
    row.start <- row.start + add
    row.stop  <- row.start + ref.length - 1

    #fill in matrix
    result[row.start:row.stop, ] <- data
  }
  return(result)
}

determine_matrix_size <- function(conditions,
                                  ref.length,
                                  total.indices,
                                  total.experiments,
                                  read.length.range) {
  # helper function that determines the size of the matrix
  # Returns:
  # Matrix of the appropriate size based on the read lengths,
  # number of experiments, number of reference transcripts, sum.transcripts
  # and aggregate
  length <- conditions[["length"]]
  transcript <- conditions[["transcript"]]
  if (length & transcript) {
    return (matrix(0L,
                   nrow = total.experiments,
                   ncol = total.indices))
  } else if (length) { #condense across lengths only
    return (matrix(0L,
                   nrow = ref.length * total.experiments,
                   ncol = total.indices))
  } else if (transcript) { #condense across transcript only
    return (matrix(nrow = read.length.range * total.experiments,
                   ncol = total.indices))
  } #!transcript & !length
  return (matrix(nrow = ref.length * read.length.range * total.experiments,
                 ncol = total.indices))
}

make_datatable <- function(ref.names,
                           transcript,
                           length,
                           matched.list,
                           range.lower,
                           range.upper,
                           matrix) {
  # helper function that creates a polished data table out of the matrix
  #
  # Returns:
  # Data table containing the metagene data based on the specifications of the
  # experiments to include, read length ranges, sum.transcript, aggregate,
  # and the matrix structure

  total.list <- length(matched.list)
  ref.length <- length(ref.names)
  num.reads  <- range.upper - range.lower + 1

  #determine columns for each separate case
  if (transcript & length) {
    experiment.list   <- matched.list
    return (data.table(experiment = matched.list,
                       matrix))
  } else if (transcript) { #sum.transcripts only
    experiment.column <- rep(matched.list, each = num.reads)
    read.column       <- rep(c(range.lower:range.upper), total.list)
    return (data.table(experiment = experiment.column,
                       length     = read.column,
                       matrix))
  } else if (length) { #aggregate only
    experiment.column <- rep(matched.list, each = ref.length)
    transcript.column <- rep(ref.names, total.list)
    return (data.table(experiment = experiment.column,
                       transcript = transcript.column,
                       matrix))
  }
  #!sum.transcripts and !aggregate
  experiment.column <- rep(matched.list, each = num.reads * ref.length)
  transcripts       <- rep(ref.names, total.list * num.reads)
  ref.read          <- rep(c(range.lower:range.upper), each = ref.length)
  read.length       <- rep(ref.read, total.list)
  return (data.table(experiment  = experiment.column,
                     transcript  = transcripts,
                     read.length = read.length,
                     matrix))
}
