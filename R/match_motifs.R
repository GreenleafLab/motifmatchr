
#' motifMatches
#'
#' get motif matches from SummarizedExperiment object
#' @param object SummarizedExperiment object with matches assay
#' @return matrix with scores
#' @export
#' @examples
#'
#' data(example_motifs, package = "motifmatchr")
#'
#' # Make a set of peaks
#' peaks <- GenomicRanges::GRanges(seqnames = c("chr1","chr2","chr2"),
#'                 ranges = IRanges::IRanges(start = c(76585873,42772928,
#'                                           100183786),
#'                                           width = 500))
#'
#' # Get motif matches for example motifs
#' motif_ix <- matchMotifs(example_motifs, peaks,
#'                         genome = "BSgenome.Hsapiens.UCSC.hg19")
#'
#' motifMatches(motif_ix)
setGeneric("motifMatches", function(object) standardGeneric("motifMatches"))

#' @describeIn motifMatches method for SummarizedExperiment
#' @export
setMethod("motifMatches", c(object = "SummarizedExperiment"),
          function(object) {
              if ("motifMatches" %ni% assayNames(object)){
                  if ("motif_matches" %in% assayNames(object)){
                      warning("motif_matches assay name deprecated",
                              "motifmatchr now using motifMatches")
                      return(assays(object)$motif_matches)
                  }
                  stop("No motifMatches in this object")
              }
              assays(object)$motifMatches
          })


#' motifScores
#'
#' get motif scores from SummarizedExperiment object
#' @param object SummarizedExperiment object with scores assay
#' @return matrix with scores
#' @export
#' @examples
#'
#' data(example_motifs, package = "motifmatchr")
#'
#' # Make a set of peaks
#' peaks <- GenomicRanges::GRanges(seqnames = c("chr1","chr2","chr2"),
#'                 ranges = IRanges::IRanges(start = c(76585873,42772928,
#'                                           100183786),
#'                                           width = 500))
#'
#' # Get motif matches for example motifs
#' motif_ix <- matchMotifs(example_motifs, peaks,
#'                          genome = "BSgenome.Hsapiens.UCSC.hg19",
#'                          out = "scores")
#'
#' motifScores(motif_ix)
setGeneric("motifScores", function(object) standardGeneric("motifScores"))

#' @describeIn motifScores method for SummarizedExperiment
#' @export
setMethod("motifScores", c(object = "SummarizedExperiment"), function(object) {
    if ("motifScores" %ni% assayNames(object)){
        if ("motif_scores" %in% assayNames(object)){
            warning("motif_scores assay name deprecated",
                    "motifmatchr now using motifScores")
            return(assays(object)$motif_scores)
        }
        stop("No motifScores in this object")
    }
    assays(object)$motifScores
})


#' motifCounts
#'
#' get motif counts from SummarizedExperiment object
#' @param object SummarizedExperiment object with counts assay
#' @return matrix with counts
#' @export
#' @examples
#'
#' data(example_motifs, package = "motifmatchr")
#'
#' # Make a set of peaks
#' peaks <- GenomicRanges::GRanges(seqnames = c("chr1","chr2","chr2"),
#'                 ranges = IRanges::IRanges(start = c(76585873,42772928,
#'                                           100183786),
#'                                           width = 500))
#'
#' # Get motif matches for example motifs
#' motif_ix <- matchMotifs(example_motifs, peaks,
#'                          genome = "BSgenome.Hsapiens.UCSC.hg19",
#'                          out = "scores")
#'
#' motifCounts(motif_ix)
setGeneric("motifCounts", function(object) standardGeneric("motifCounts"))

#' @describeIn motifCounts method for SummarizedExperiment
#' @export
setMethod("motifCounts", c(object = "SummarizedExperiment"), function(object) {
    if ("motifCounts" %ni% assayNames(object)){
        if ("motif_counts" %in% assayNames(object)){
            warning("motif_counts assay name deprecated",
                    "motifmatchr now using motifCounts")
            return(assays(object)$motif_counts)
        }
        stop("No motifCounts in this object")
    }
    assays(object)$motifCounts
})

matchMotifs_helper <- function(pwms, seqs, bg, p.cutoff, w, out, ranges) {

    motif_mats <- convert_pwms(pwms, bg)

    if (out == "matches") {
        tmp_out <- get_motif_ix(motif_mats, seqs, bg, p.cutoff, w)
        if (is.null(ranges)) {
            out <- SummarizedExperiment(assays =
                                            list(motifMatches = as(tmp_out,
                                                              "lMatrix")),
                                        colData =
                                            DataFrame(name = name(pwms),
                                                      row.names = names(pwms)))
        } else {
            out <- SummarizedExperiment(assays =
                                            list(motifMatches = as(tmp_out,
                                                              "lMatrix")),
                                        rowRanges = ranges,
                                        colData =
                                            DataFrame(name = name(pwms),
                                                      row.names = names(pwms)))
        }
    } else if (out == "scores") {
        tmp_out <- get_motif_ix_plus(motif_mats, seqs, bg, p.cutoff, w)
        tmp_out$motifMatches <- as(tmp_out$motifMatches, "lMatrix")
        if (is.null(ranges)) {
            out <- SummarizedExperiment(assays = tmp_out,
                                        colData =
                                            DataFrame(name = name(pwms),
                                                      row.names = names(pwms)))
        } else {
            out <- SummarizedExperiment(assays = tmp_out,
                                        rowRanges = ranges,
                                        colData =
                                            DataFrame(name = name(pwms),
                                                      row.names = names(pwms)))
        }
    } else {
        tmp_out <- get_motif_positions(motif_mats, seqs, bg, p.cutoff, w)
        if (is.null(ranges)) {
            out <- lapply(1:length(motif_mats), function(x) {
                m_ix <- which(tmp_out$motif_ix == x - 1)
                IRangesList(lapply(seq_along(seqs),
                       function(y){
                           ms_ix <- m_ix[which(tmp_out$seq_ix[m_ix] == y -1)]
                           tmp <- IRanges(start = tmp_out$pos[ms_ix] + 1,
                               width = ncol(motif_mats[[x]]))
                            mcols(tmp) <- DataFrame(strand = tmp_out$strand[ms_ix],
                                        score = tmp_out$score[ms_ix])
                            tmp
                       }))
            })
            names(out) <- names(pwms)
        } else {
            out <- lapply(1:length(motif_mats), function(x) {
                m_ix <- which(tmp_out$motif_ix == x - 1)
                GRanges(seqnames(ranges)[tmp_out$seq_ix[m_ix] +  1],
                        IRanges(start =
                                    start(ranges[tmp_out$seq_ix[m_ix] + 1]) +
                                    tmp_out$pos[m_ix],
                                width = ncol(motif_mats[[x]])),
                        strand = tmp_out$strand[m_ix],
                        score = tmp_out$score[m_ix])
            })
            names(out) <- names(pwms)
            out <- GRangesList(out)
        }
    }
    return(out)
}



#' matchMotifs
#'
#' Find motif matches
#' @param pwms either \code{\link[TFBSTools]{PFMatrix}},
#' \code{\link[TFBSTools]{PFMatrixList}}, \code{\link[TFBSTools]{PWMatrix}},
#' \code{\link[TFBSTools]{PWMatrixList}}
#' @param subject either \code{\link[GenomicRanges]{GenomicRanges}},
#' \code{\link[Biostrings]{DNAStringSet}}, \code{\link[Biostrings]{DNAString}},
#' or character vector
#' @param genome BSgenome object, \code{\link[Biostrings]{DNAStringSet}}, or
#' \code{\link[Rsamtools]{FaFile}}, or short string signifying genome build
#' recognized by \code{\link[BSgenome]{getBSgenome}}. Only required if
#' subject is \code{\link[GenomicRanges]{GenomicRanges}} or
#' \code{\link[SummarizedExperiment]{RangedSummarizedExperiment}} or if bg is set
#' to "genome"
#' @param bg background nucleotide frequencies. Default is to compute based on
#' subject, i.e. the specific set of sequences being evaluated. See Details.
#' @param out what to return? see return section
#' @param p.cutoff p-value cutoff for returning motifs
#' @param w parameter controlling size of window for filtration; default is 7
#' @param ranges if subject is not GenomicRanges or RangedSummarizedExperiment,
#'  these ranges can be used to specify what ranges the input sequences
#'  correspond to. These ranges will be incorporated into the
#'  SummarizedExperiment output if out is "matches" or "scores" or will be used
#'  to give absolute positions of motifs if out is "positions"
#' @param ... additional arguments depending on inputs
#' @details Background nucleotide frequencies can be set to "subject" to use the
#' subject sequences or ranges for computing the nucleotide frequencies,
#' "genome" for using the genomice frequencies (in which case a genome must be
#' specified), "even" for using 0.25 for each base, or a numeric vector with A,
#' C, G, and T frequencies.
#' @return Either returns a SummarizedExperiment with a sparse matrix with
#'  values set to TRUE for a match (if out == 'matches'), a
#'  SummarizedExperiment with a matches matrix as well as matrices with the
#'  maximum motif score and total motif counts (if out == 'scores'), or a
#'  \code{\link[GenomicRanges]{GenomicRangesList}} or a list of
#'  \code{\link[IRanges]{IRangesList}} with all the positions of matches
#'  (if out == 'positions')
#' @export
#' @examples
#'
#' data(example_motifs, package = "motifmatchr")
#'
#' # Make a set of peaks
#' peaks <- GenomicRanges::GRanges(seqnames = c("chr1","chr2","chr2"),
#'                 ranges = IRanges::IRanges(start = c(76585873,42772928,
#'                                           100183786),
#'                                           width = 500))
#'
#' # Get motif matches for example motifs
#' motif_ix <- matchMotifs(example_motifs, peaks, genome = "BSgenome.Hsapiens.UCSC.hg19")
#'
setGeneric("matchMotifs",
           function(pwms, subject, ...) standardGeneric("matchMotifs"))

#' @describeIn matchMotifs PWMatrixList/DNAStringSet
#' @export
setMethod("matchMotifs", signature(pwms = "PWMatrixList",
                                    subject = "DNAStringSet"),
          function(pwms,
                   subject,
                   genome = NULL,
                   bg = c("subject","genome","even"),
                   out = c("matches", "scores", "positions"),
                   p.cutoff = 5e-05, w = 7, ranges = NULL) {
              out <- match.arg(out)

              if (is.numeric(bg)){
                  bg <- check_bg(bg)
              } else{
                  bg_method <- match.arg(bg)
                  bg <- get_bg(bg_method, subject, genome)
              }

              seqs <- as.character(subject)

              matchMotifs_helper(pwms, seqs, bg, p.cutoff, w, out, ranges)
          })

#' @describeIn matchMotifs PWMatrixList/character
#' @export
setMethod("matchMotifs", signature(pwms = "PWMatrixList",
                                    subject = "character"),
          function(pwms,
                   subject,
                   genome = NULL,
                   bg = c("subject","genome","even"),
                   out = c("matches", "scores",  "positions"),
                   p.cutoff = 5e-05, w = 7, ranges = NULL) {

              out <- match.arg(out)

              if (is.numeric(bg)){
                  bg <- check_bg(bg)
              } else{
                  bg_method <- match.arg(bg)
                  bg <- get_bg(bg_method, subject, genome)
              }

              matchMotifs_helper(pwms, subject, bg, p.cutoff, w, out, ranges)
          })

#' @describeIn matchMotifs PWMatrixList/DNAString
#' @export
setMethod("matchMotifs", signature(pwms = "PWMatrixList",
                                    subject = "DNAString"),
          function(pwms,
                   subject,
                   genome = NULL,
                   bg = c("subject","genome","even"),
                   out = c("matches", "scores", "positions"),
                   p.cutoff = 5e-05, w = 7, ranges = NULL) {

              out <- match.arg(out)

              if (is.numeric(bg)){
                  bg <- check_bg(bg)
              } else{
                  bg_method <- match.arg(bg)
                  bg <- get_bg(bg_method, subject, genome)
              }

              seqs <- as.character(subject)

              matchMotifs_helper(pwms, seqs, bg, p.cutoff, w, out, ranges)
          })

#' @describeIn matchMotifs PWMatrixList/GenomicRanges
#' @export
setMethod("matchMotifs", signature(pwms = "PWMatrixList",
                                    subject = "GenomicRanges"),
          function(pwms,
                   subject,
                   genome = GenomeInfoDb::genome(subject),
                   bg = c("subject","genome","even"),
                   out = c("matches", "scores", "positions"),
                   p.cutoff = 5e-05, w = 7) {
              out <- match.arg(out)
              GenomicRanges::strand(subject) <- "+"
              genome <- validate_genome_input(genome)
              seqs <- getSeq(genome, subject)

              if (is.numeric(bg)){
                  bg <- check_bg(bg)
              } else{
                  bg_method <- match.arg(bg)
                  bg <- get_bg(bg_method, seqs, genome)
              }

              seqs <- as.character(seqs)

              matchMotifs_helper(pwms, seqs, bg, p.cutoff, w, out, subject)
          })

#' @describeIn matchMotifs PWMatrixList/RangedSummarizedExperiment
#' @export
setMethod("matchMotifs", signature(pwms = "PWMatrixList",
                                    subject = "RangedSummarizedExperiment"),
          function(pwms, subject,
                   genome = GenomeInfoDb::genome(subject),
                   bg = c("subject","genome","even"),
                   out = c("matches", "scores", "positions"),
                   p.cutoff = 5e-05, w = 7) {
              out <- match.arg(out)
              matchMotifs(pwms, rowRanges(subject),
                           genome = genome,
                           bg = bg,
                           out = out,
                           p.cutoff = p.cutoff,
                           w = w)
          })

#' @describeIn matchMotifs PWMatrixList/BSGenomeViews
#' @export
setMethod("matchMotifs", signature(pwms = "PWMatrixList",
                                    subject = "BSgenomeViews"),
          function(pwms,
                   subject,
                   bg = c("subject","genome","even"),
                   out = c("matches", "scores", "positions"),
                   p.cutoff = 5e-05, w = 7) {

              out <- match.arg(out)
              seqs <- as.character(subject)
              ranges <- BSgenome::granges(subject)

              if (is.numeric(bg)){
                  bg <- check_bg(bg)
              } else{
                  bg_method <- match.arg(bg)
                  bg <- get_bg(bg_method, subject, BSgenome::subject(subject))
              }

              seqs <- as.character(subject)

              matchMotifs_helper(pwms, seqs, bg, p.cutoff, w, out, ranges)
          })

### PFMatrixList ---------------------------------------------------------------


#' @describeIn matchMotifs PFMatrixList/ANY
#' @export
setMethod("matchMotifs", signature(pwms = "PFMatrixList", subject = "ANY"),
          function(pwms,
                   subject,
                   ...) {


              pwms_list <- do.call(PWMatrixList, lapply(pwms, toPWM))
              matchMotifs(pwms_list,
                           subject,
                           ...)
          })

# Single PWM input -------------------------------------------------------------

#' @describeIn matchMotifs PWMatrix/ANY
#' @export
setMethod("matchMotifs", signature(pwms = "PWMatrix", subject = "ANY"),
          function(pwms,
                   subject,
                   ...) {

              pwms_list <- PWMatrixList(pwms)
              matchMotifs(pwms_list,
                           subject,
                           ...)
          })


# Single PFM -------------------------------------------------------------------

#' @describeIn matchMotifs PFMatrix/ANY
#' @export
setMethod("matchMotifs",
          signature(pwms = "PFMatrix", subject = "ANY"),
          function(pwms,
                   subject,
                   ...) {

              pwms_list <- PWMatrixList(toPWM(pwms))
              matchMotifs(pwms_list,
                           subject,
                           ...)
          })
