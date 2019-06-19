#' Information on the RNA-Seq data of the experiments, if any
#'
#' \code{\link{get_rnaseq}} returns a data.table containing information on the transcript name, experiment, and
#' sequence abundance
#'
#' As a default value, experiment.list is presumed to include all of the
#' experiments within a ribo file. RNA-Seq data is an optional dataset to
#' include in a .ribo file. The experiments in experiment.list are checked
#' for experiment existence in the ribo file and then checked for RNA-seq data.
#'
#' The returned data.table can either be in the tidy format for easier data
#' cleaning or in a condensed non-tidy format. The data will present RNA-seq abundance
#' for each transcript in each valid experiment in experiment.list.
#' @examples
#' #generate the ribo object
#' file.path <- system.file("extdata", "sample.ribo", package = "ribor")
#' sample <- ribo(file.path)
#'
#' #list out the experiments of interest that have RNA-Seq data
#' experiments = c("Hela_1", "Hela_2", "WT_1")
#' rnaseq.data <- get_rnaseq(ribo.object = sample,
#'                           tidy = TRUE,
#'                           experiments = experiments)
#'
#' @param ribo.object S3 object of class "ribo"
#' @param tidy logical value denoting whether or not the user wants a tidy format
#' @param experiments list of experiment names
#'
#' @return
#' Returns a data table that contains the transcript name, experiment, and
#' RNA-seq abundance
#'
#' @seealso \code{\link{ribo}} to generate the necessary ribo.object parameter
#'
#' @importFrom rhdf5 h5ls h5read
#' @importFrom data.table data.table
#' @importFrom tidyr gather
#' @export
get_rnaseq <- function(ribo.object,
                       tidy = TRUE,
                       experiments = get_experiments(ribo.object)) {

  ribo.experiments <- get_experiments(ribo.object)


  check_experiments(ribo.object, experiments)

  #get just the experiments that exist
  rnaseq.experiments <- check_rnaseq(ribo.object, experiments)

  #generate appropriate matrix size in the untidy version
  ref.names <- get_reference_names(ribo.object)
  ref.length <- length(ref.names)
  total.experiments <- length(rnaseq.experiments)
  num.regions <- 5
  handle <- ribo.object$handle

  result <- matrix(nrow = ref.length * total.experiments,
                   ncol = num.regions)
  colnames(result) <- c("UTR5", "UTR5J", "CDS", "UTR3J", "UTR3")

  #generate the untidy version
  for (index in 1:total.experiments) {
      experiment <- rnaseq.experiments[index]
      path = paste("/experiments/", experiment, "/rnaseq/rnaseq", sep = "")
      row.start <- (index - 1) * ref.length + 1
      row.stop <- row.start + ref.length - 1
      result[row.start:row.stop, ] <- t(h5read(handle, path))
  }

  rnaseq <- rep(rnaseq.experiments, each = ref.length)
  transcripts <- rep(ref.names, total.experiments)

  result <- data.table(experiment = rnaseq,
                       transcript = transcripts,
                       result)
  if (tidy) {
    result <- gather(result, "region", "abundance", c("UTR5":"UTR3"))
  }
  return(result)
}

check_rnaseq <- function(ribo.object, experiments) {
  #obtain the rnaseq data
  rnaseq <- get_info(ribo.object)[["contents"]][, c("names", "rna.seq")]
  has.rnaseq <- rnaseq[rnaseq$rna.seq == TRUE, ]
  has.rnaseq <- has.rnaseq$names
  check_experiments(ribo.object, experiments)

  #find the experiments in the experiment.list that do not have coverage and print warnings
  check <- setdiff(experiments, has.rnaseq)
  if (length(check)) {
    for (experiment in check) {
      warning("'", experiment, "'", " did not have RNA-Seq data.", call. = FALSE)
    }
    warning("Param experiments contains experiments that did not have RNA-Seq data.
            The return value ignores these experiments.", call. = FALSE)
  }
  return(intersect(experiments, has.rnaseq))
}
