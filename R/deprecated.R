#' Deprecated functions in motifmatchr
#'
#' motifmatchr has moved functions and methods to camelCase from snake_case.
#' The following functions have been deprecated and replaced with a different
#' name:
#' \itemize{
#'   \item motif_matches is now \code{\link{motifMatches}}
#'   \item motif_counts is now \code{\link{motifCounts}}
#'   \item motif_scores is now  \code{\link{motifScores}}
#'   \item match_motifs is now  \code{\link{matchMotifs}}
#' }
#' @author Alicia Schep
#' @rdname motifmatchr_deprecated
#' @param ... arguments passed to new function
#' @return calls the replacement function
#' @name motifmatchr_deprecated
NULL

#' @rdname motifmatchr_deprecated
#' @export
motif_matches <- function(...)
{
    .Deprecated("motifMatches")
    motifMatches(...)
}

#' @rdname motifmatchr_deprecated
#' @export
motif_counts <- function(...)
{
    .Deprecated("motifCounts")
    motifCounts(...)
}

#' @rdname motifmatchr_deprecated
#' @export
motif_scores <- function(...)
{
    .Deprecated("motifScores")
    motifScores(...)
}

#' @rdname motifmatchr_deprecated
#' @export
match_motifs <- function(...)
{
    .Deprecated("matchMotifs")
    matchMotifs(...)
}

