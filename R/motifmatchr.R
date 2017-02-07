#' motifmatchr: Fast Motif Matching in R
#'
#' The motifmatchr package is designed for analyzing many sequences and many
#' motifs to find which sequences contain which motifs.
#'
#' motifmatchr uses the MOODS C++ library (developedby Pasi Rastas,
#' Janne Korhonen, and Petri Martinmaki) internally for motif matching.
#'
#' The primary method of motifmatchr is \code{\link{match_motifs}}, which takes
#' in motif PWMs/PFMs and genomic ranges or sequences and returns either which
#' ranges/sequences match which motifs or the positions of the matches.
#'
#' Compared with alternative motif matching functions available in Bioconductor
#' (e.g. matchPWM in Biostrings or searchSeq in TFBSTools), motifmatchr is
#' designed specifically for the use case of determining whether many different
#' sequences/ranges contain many different motifs.
#'
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
#' @author Alicia Schep
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
