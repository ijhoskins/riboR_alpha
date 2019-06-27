get_reference_names <- function(ribo.object){
  # Retrieves the reference transcript names 
  check_ribo(ribo.object)
  return(h5read(ribo.object$handle&'reference', name = "reference_names"))
}


get_reference_lengths <- function(ribo.object){
  # Retrieves the reference transcript lengths 
  check_ribo(ribo.object)
  row.names <- h5read(ribo.object$handle&'reference', 
                      name = "reference_names")
  lengths   <- h5read(ribo.object$handle&'reference', 
                      name = "reference_lengths")
  return(data.table(transcript = row.names, length = lengths))
}


get_content_info <- function(ribo.handle) {
  experiments        <- h5ls(ribo.handle&'experiments', recursive = FALSE)$name
  length <- length(experiments)
  
  #creates the separate lists for reads, coverage, rna.seq, and metadata
  #to eventually put in a data frame
  reads.list    <- vector(mode = "integer", length = length)
  coverage.list <- vector(mode = "logical", length = length)
  rna.seq.list  <- vector(mode = "logical", length = length)
  metadata.list <- vector(mode = "logical", length = length)
  
  #ls function provides information about the contents of each experiment
  ls <- h5ls(ribo.handle)
  
  #loop over all of the experiments
  for (i in 1:length) {
    experiment <- experiments[i]
    #gathers information on the number of reads for each experiment by looking at
    #the attributes
    name           <- paste("/experiments/", experiment, sep = "")
    attribute      <- h5readAttributes(ribo.handle, name)
    reads.list[i]     <- attribute[["total_reads"]]
    
    #creates separate logical lists to denote the presence of
    #reads, coverage, RNA-seq, metadata
    metadata.list[i]  <- ("metadata" %in% names(attribute))
    
    group.contents <- ls[ls$group == name,]
    group.names    <- group.contents$name
    
    coverage.list[i]  <- ("coverage" %in% group.names)
    rna.seq.list[i]   <- ("rnaseq" %in% group.names)
  }
  
  experiments.info       <- data.table(experiment  = experiments,
                                       total.reads = reads.list,
                                       coverage    = coverage.list,
                                       rna.seq     = rna.seq.list,
                                       metadata    = metadata.list)
  return(experiments.info)
}


get_attributes <- function(ribo.object) {
  # Retrieves the attributes of the ribo.object 
  handle  <- ribo.object$handle
  attribute <- h5readAttributes(handle, "/")
  return(attribute[-which(names(attribute) == "time")])
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
