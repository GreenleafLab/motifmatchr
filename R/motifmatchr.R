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

#' example_motifs
#'
#' A few example motifs from JASPAR 2016 for trying out motifmatchr
#' @docType data
#' @keywords datasets
#' @name example_motifs
#' @usage data(example_motifs)
#' @format \code{\link[TFBSTools]{PFMatrixList}} of length 3
NULL



.onUnload <- function(libpath) {
  library.dynam.unload("motifmatchr", libpath)
}
