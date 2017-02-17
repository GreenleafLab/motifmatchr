# Not exported

# Not-in operator --------------------------------------------------------------

"%ni%" <- Negate("%in%")

## get_nuc_freqs function to find background sequence --------------------------

setGeneric("get_nuc_freqs", function(subject, ...) standardGeneric("get_nuc_freqs"))

setMethod("get_nuc_freqs", signature(subject = "DNAStringSet"),
          function(subject) {
            nucFreqs <- colSums(letterFrequency(subject, c("A", "C", "G", "T")))
            nucFreqs <- nucFreqs/sum(nucFreqs)
            return(nucFreqs)
          })

setMethod("get_nuc_freqs", signature(subject = "DNAString"),
          function(subject) {
            nucFreqs <- letterFrequency(subject, c("A", "C", "G", "T"))
            nucFreqs <- nucFreqs/sum(nucFreqs)
            return(nucFreqs)
          })

setMethod("get_nuc_freqs", signature(subject = "character"),
          function(subject) {
            if (length(subject) == 1) {
              return(get_nuc_freqs(DNAString(subject)))
            } else {
              return(get_nuc_freqs(DNAStringSet(subject)))
            }
          })


## convert_pwm to adjust background in pwm model ---------------------------------

convert_pwms <- function(pwms, bg_freqs) {
  stopifnot(inherits(pwms, "PWMatrixList"))
  lapply(pwms, convert_pwm, bg_freqs)
}



check_bg <- function(bg_freqs){
  if (length(bg_freqs) != 4) stop("Invalid background frequencies -- should be length 4")
  if (!all.equal(sum(bg_freqs),1) || min(bg_freqs) <= 0){
    stop("Invalid background frequencies. Should sum to 1")
  }

  if (!is.null(names(bg_freqs))){
    if (!all(names(c("A","C","G","T") %in% bg_freqs))){
      stop("Background nucleotide frequencies have names that don't match nucleotides! (A,C,G,T)")
    } else{
      bg_freqs <- bg_freqs[c("A","C","G","T")]
    }
  }
  return(bg_freqs)
}

convert_pwm <- function(pwm, bg_freqs) {
  type <- pwm_type(pwm)
  out <- as.matrix(pwm)
  if (type == "prob") {
    norm_mat <- matrix(bg_freqs, nrow = 4, ncol = length(pwm), byrow = FALSE)
    out <- log(as.matrix(pwm)/norm_mat)
  } else if (type == "log2") {
    norm_mat <- matrix(log2(bg(pwm)) - log2(bg_freqs), nrow = 4,
                       ncol = length(pwm),
      byrow = FALSE)
    out <- as.matrix(pwm) - norm_mat
  } else if (type == "log") {
    norm_mat <- matrix(log(bg(pwm)) - log(bg_freqs), nrow = 4,
                       ncol = length(pwm),
      byrow = FALSE)
    out <- as.matrix(pwm) - norm_mat
  }
  return(out)
}

#' pwm_type
#'
#' Determines type of PWM
#' @param pwm PWMatrix object
#' @return 'log','log2', or 'frequency' depending on type of pwm
#' @export
#' @keywords internal
#' @examples
#'
#' data(example_motifs, package = "motifmatchr")
#' pwm_type(TFBSTools::toPWM(example_motifs[[1]]))
#' pwm_type(TFBSTools::toPWM(example_motifs[[1]], type = "prob"))
pwm_type <- function(pwm) {
  # Determine whether un-logged, natural log, or log2
  if (isTRUE(all.equal(colSums(as.matrix(pwm)), rep(1, length(pwm))))) {
    return("frequency")
  } else if (isTRUE(all.equal(colSums(2^(as.matrix(pwm)) *
                                      matrix(bg(pwm),
                                             byrow = FALSE,
                                             ncol = length(pwm),
                                             nrow = 4)),
                              rep(1, length(pwm)), tolerance = 10^-5))) {
    return("log2")
  } else if (isTRUE(all.equal(colSums(exp(as.matrix(pwm)) *
                                      matrix(bg(pwm),
                                             byrow = FALSE,
                                             ncol = length(pwm),
                                             nrow = 4)),
                              rep(1, length(pwm)), tolerance = 10^-5))) {
    return("log")
  } else {
    stop("Can't determine format of PWM -- should be numeric frequency summing to 1 or log or log2 odds ratio")
  }
}

