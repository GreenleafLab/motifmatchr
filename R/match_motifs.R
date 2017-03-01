#' motif_matches
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
#' # Get motif matches for example motifs  (using hg19 genome, the default)
#' motif_ix <- match_motifs(example_motifs, peaks)
#'
#' motif_matches(motif_ix)
setGeneric("motif_matches", function(object) standardGeneric("motif_matches"))

#' @describeIn motif_matches method for SummarizedExperiment
#' @export
setMethod("motif_matches", c(object = "SummarizedExperiment"),
          function(object) {
              if ("matches" %ni% assayNames(object))
                  stop("No matches in this object")
              assays(object)$matches
          })


#' motif_scores
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
#' # Get motif matches for example motifs  (using hg19 genome, the default)
#' motif_ix <- match_motifs(example_motifs, peaks, out = "scores")
#'
#' motif_scores(motif_ix)
setGeneric("motif_scores", function(object) standardGeneric("motif_scores"))

#' @describeIn motif_scores method for SummarizedExperiment
#' @export
setMethod("motif_scores", c(object = "SummarizedExperiment"), function(object) {
    if ("scores" %ni% assayNames(object))
        stop("No scores in this object")
    assays(object)$scores
})


#' motif_counts
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
#' # Get motif matches for example motifs  (using hg19 genome, the default)
#' motif_ix <- match_motifs(example_motifs, peaks, out = "scores")
#'
#' motif_counts(motif_ix)
setGeneric("motif_counts", function(object) standardGeneric("motif_counts"))

#' @describeIn motif_counts method for SummarizedExperiment
#' @export
setMethod("motif_counts", c(object = "SummarizedExperiment"), function(object) {
    if ("counts" %ni% assayNames(object))
        stop("No counts in this object")
    assays(object)$counts
})

match_motifs_helper <- function(pwms, seqs, bg, p.cutoff, w, out, ranges) {

    motif_mats <- convert_pwms(pwms, bg)

    if (out == "matches") {
        tmp_out <- get_motif_ix(motif_mats, seqs, bg, p.cutoff, w)
        if (is.null(ranges)) {
            out <- SummarizedExperiment(assays =
                                            list(matches = as(tmp_out,
                                                              "lMatrix")),
                                        colData =
                                            DataFrame(pwm = pwms,
                                                      name = name(pwms),
                                                      row.names = names(pwms)))
        } else {
            out <- SummarizedExperiment(assays =
                                            list(matches = as(tmp_out,
                                                              "lMatrix")),
                                        rowRanges = ranges,
                                        colData =
                                            DataFrame(pwm = pwms,
                                                      name = name(pwms),
                                                      row.names = names(pwms)))
        }
    } else if (out == "scores") {
        tmp_out <- get_motif_ix_plus(motif_mats, seqs, bg, p.cutoff, w)
        tmp_out$matches <- as(tmp_out$matches, "lMatrix")
        if (is.null(ranges)) {
            out <- SummarizedExperiment(assays = tmp_out,
                                        colData =
                                            DataFrame(pwm = pwms,
                                                      name = name(pwms),
                                                      row.names = names(pwms)))
        } else {
            out <- SummarizedExperiment(assays = tmp_out,
                                        rowRanges = ranges,
                                        colData =
                                            DataFrame(pwm = pwms,
                                                      name = name(pwms),
                                                      row.names = names(pwms)))
        }
    } else {
        tmp_out <- get_motif_positions(motif_mats, seqs, bg, p.cutoff, w)
        if (is.null(ranges)) {
            out <- lapply(1:length(motif_mats), function(x) {
                m_ix <- which(tmp_out$motif_ix == x - 1)
                tmp <- IRanges(start = tmp_out$pos[m_ix] + 1,
                               width = ncol(motif_mats[[x]]))
                mcols(tmp) <- DataFrame(strand = tmp_out$strand[m_ix],
                                        score = tmp_out$score[m_ix])
                tmp
            })
            names(out) <- names(motif_mats)
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
        }
        names(out) <- names(pwms)
    }
    return(out)
}


#' match_motifs
#'
#' Find motif matches
#' @param pwms either \code{\link[TFBSTools]{PFMatrix}},
#' \code{\link[TFBSTools]{PFMatrixList}}, \code{\link[TFBSTools]{PWMatrix}},
#' \code{\link[TFBSTools]{PWMatrixList}}
#' @param subject either \code{\link[GenomicRanges]{GenomicRanges}},
#' \code{\link[Biostrings]{DNAStringSet}}, \code{\link[Biostrings]{DNAString}},
#' or character vector
#' @param genome BSgenome object, only used if subect is
#' \code{\link[GenomicRanges]{GenomicRanges}}
#' @param bg background nucleotide frequencies. if not provided, computed from
#' subject
#' @param out what to return? see return section
#' @param p.cutoff p-value cutoff for returning motifs
#' @param w parameter controlling size of window for filtration; default is 7
#' @param ranges if subject is not GenomicRanges, ranges to use when out is
#' positions
#' @param ... additional arguments depending on inputs
#' @return Either returns a SummarizedExperiment with a sparse matrix with
#'  values set to TRUE for a match (if return == 'matches'), a
#'  SummarizedExperiment with a matches matrix as well as matrices with the
#'  maximum motif score and total motif counts (if return == 'scores'), or a
#'  GenomicRanges or IRanges object with all the positions of matches or
#'  \code{\link[GenomicRanges]{GenomicRanges}} if positions
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
#' # Get motif matches for example motifs  (using hg19 genome, the default)
#' motif_ix <- match_motifs(example_motifs, peaks)
#'
setGeneric("match_motifs",
           function(pwms, subject, ...) standardGeneric("match_motifs"))

#' @describeIn match_motifs PWMatrixList/DNAStringSet
#' @export
setMethod("match_motifs", signature(pwms = "PWMatrixList",
                                    subject = "DNAStringSet"),
          function(pwms, subject, genome = NULL, bg = NULL,
                   out = c("matches", "scores", "positions"),
                   p.cutoff = 5e-05, w = 7, ranges = NULL) {
              out <- match.arg(out)

              if (is.null(bg)) {
                  bg <- get_nuc_freqs(subject)
              } else {
                  bg <- check_bg(bg)
              }
              seqs <- as.character(subject)

              match_motifs_helper(pwms, seqs, bg, p.cutoff, w, out, ranges)
          })

#' @describeIn match_motifs PWMatrixList/character
#' @export
setMethod("match_motifs", signature(pwms = "PWMatrixList",
                                    subject = "character"),
          function(pwms, subject, genome = NULL, bg = NULL,
                   out = c("matches", "scores",  "positions"),
                   p.cutoff = 5e-05, w = 7, ranges = NULL) {

              out <- match.arg(out)

              if (is.null(bg)) {
                  bg <- get_nuc_freqs(subject)
              } else {
                  bg <- check_bg(bg)
              }
              match_motifs_helper(pwms, subject, bg, p.cutoff, w, out, ranges)
          })

#' @describeIn match_motifs PWMatrixList/DNAString
#' @export
setMethod("match_motifs", signature(pwms = "PWMatrixList",
                                    subject = "DNAString"),
          function(pwms, subject, genome = NULL, bg = NULL,
                   out = c("matches", "scores", "positions"),
                   p.cutoff = 5e-05, w = 7, ranges = NULL) {
              out <- match.arg(out)
              if (is.null(bg)) {
                  bg <- get_nuc_freqs(subject)
              } else {
                  bg <- check_bg(bg)
              }

              seqs <- as.character(subject)

              match_motifs_helper(pwms, seqs, bg, p.cutoff, w, out, ranges)
          })

#' @describeIn match_motifs PWMatrixList/GenomicRanges
#' @export
setMethod("match_motifs", signature(pwms = "PWMatrixList",
                                    subject = "GenomicRanges"),
          function(pwms, subject, genome = BSgenome.Hsapiens.UCSC.hg19,
                   bg = NULL, out = c("matches", "scores", "positions"),
                   p.cutoff = 5e-05, w = 7,
                   ranges = NULL) {
              out <- match.arg(out)
              GenomicRanges::strand(subject) <- "+"
              seqs <- getSeq(genome, subject)
              if (is.null(bg)) {
                  bg <- get_nuc_freqs(seqs)
              } else {
                  bg <- check_bg(bg)
              }

              seqs <- as.character(seqs)

              match_motifs_helper(pwms, seqs, bg, p.cutoff, w, out, subject)
          })

#' @describeIn match_motifs PWMatrixList/RangedSummarizedExperiment
#' @export
setMethod("match_motifs", signature(pwms = "PWMatrixList",
                                    subject = "RangedSummarizedExperiment"),
          function(pwms, subject, genome = BSgenome.Hsapiens.UCSC.hg19,
                   bg = NULL, out = c("matches", "scores", "positions"),
                   p.cutoff = 5e-05, w = 7,
                   ranges = NULL) {
              out <- match.arg(out)
              match_motifs(pwms, rowRanges(subject), genome, bg, out, p.cutoff,
                           w)
          })

### PFMatrixList ---------------------------------------------------------------


#' @describeIn match_motifs PFMatrixList/ANY
#' @export
setMethod("match_motifs", signature(pwms = "PFMatrixList", subject = "ANY"),
          function(pwms,
                   subject, genome = BSgenome.Hsapiens.UCSC.hg19,
                   bg = NULL,
                   out = c("matches", "scores", "positions"),
                   p.cutoff = 5e-05, w = 7, ranges = NULL) {
              out <- match.arg(out)
              pwms_list <- do.call(PWMatrixList, lapply(pwms, toPWM))
              match_motifs(pwms_list,
                           subject,
                           genome = genome,
                           bg = bg,
                           out = out,
                           p.cutoff = p.cutoff,
                           w = w, ranges = ranges)
          })

# Single PWM input -------------------------------------------------------------

#' @describeIn match_motifs PWMatrix/ANY
#' @export
setMethod("match_motifs", signature(pwms = "PWMatrix", subject = "ANY"),
          function(pwms,
                   subject, genome = BSgenome.Hsapiens.UCSC.hg19, bg = NULL,
                   out = c("matches", "scores", "positions"), p.cutoff = 5e-05,
                   w = 7, ranges = NULL) {
              out <- match.arg(out)
              pwms_list <- PWMatrixList(pwms)
              match_motifs(pwms_list,
                           subject,
                           genome = genome,
                           bg = bg,
                           out = out,
                           p.cutoff = p.cutoff,
                           w = w, ranges = ranges)
          })


# Single PFM -------------------------------------------------------------------

#' @describeIn match_motifs PFMatrix/ANY
#' @export
setMethod("match_motifs",
          signature(pwms = "PFMatrix", subject = "ANY"),
          function(pwms,
                   subject,
                   genome = BSgenome.Hsapiens.UCSC.hg19,
                   bg = NULL,
                   out = c("matches", "scores", "positions"),
                   p.cutoff = 5e-05,
                   w = 7,
                   ranges = NULL) {
              out <- match.arg(out)
              pwms_list <- PWMatrixList(toPWM(pwms))
              match_motifs(pwms_list,
                           subject,
                           genome = genome,
                           bg = bg,
                           out = out,
                           p.cutoff = p.cutoff,
                           w = w,
                           ranges = ranges)
          })
