#' motifmatchr
#'
#' Fast motif matching in R
#' @useDynLib motifmatchr
#' @docType package
#' @import Matrix SummarizedExperiment methods
#' @importFrom S4Vectors DataFrame
#' @importFrom GenomicRanges GRanges
#' @importFrom Biostrings getSeq letterFrequency DNAString DNAStringSet
#' @importFrom BSgenome.Hsapiens.UCSC.hg19 BSgenome.Hsapiens.UCSC.hg19
#' @importFrom TFBSTools PWMatrixList toPWM name bg
#' @importFrom IRanges IRanges
#' @importMethodsFrom GenomicRanges seqnames start
#' @importMethodsFrom TFBSTools as.matrix
#' @importClassesFrom Biostrings DNAString DNAStringSet
#' @importClassesFrom GenomicRanges GenomicRanges
#' @importClassesFrom TFBSTools PWMatrix PFMatrix PWMatrixList PFMatrixList
#' @name motifmatchr
NULL
# > NULL


.onUnload <- function(libpath) {
  library.dynam.unload("motifmatchr", libpath)
}
