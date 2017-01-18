#' motifmatchr
#'
#' Fast motif matching in R
#' @useDynLib motifmatchr
#' @docType package
#' @import Matrix SummarizedExperiment
#' @importFrom S4Vectors DataFrame
#' @name motifmatchr
NULL
# > NULL


.onUnload <- function(libpath) {
  library.dynam.unload("motifmatchr", libpath)
}
