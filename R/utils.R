# Not exported

# Not-in operator --------------------------------------------------------------

"%ni%" <- Negate("%in%")

## getNucFreqs function to find background sequence --------------------------


setGeneric("getNucFreqs",
           function(subject, ...) standardGeneric("getNucFreqs"))

#' @import BSgenome
setMethod("getNucFreqs", signature(subject = "BSgenome"),
          function(subject) {
              param <- new("BSParams", X = subject, FUN = letterFrequency)
              nucFreqs <- colSums(do.call(rbind,
                                          bsapply(param,
                                                  letters =
                                                      c("A", "C", "G", "T"))))
              nucFreqs <- nucFreqs/sum(nucFreqs)
              return(nucFreqs)
          })

#' @importFrom Rsamtools scanFaIndex scanFa
setMethod("getNucFreqs", signature(subject = "FaFile"),
          function(subject) {
              nucFreqs <- c("A" = 0, "C" = 0, "G" = 0, "T" = 0)
              chroms <- scanFaIndex(subject)
              for (i in seq_along(chroms)){
                  nucFreqs <- nucFreqs +
                      letterFrequency(scanFa(subject, chroms[i])[[1]],
                                      letters = c("A", "C", "G", "T"))
              }
              nucFreqs <- nucFreqs/sum(nucFreqs)
              return(nucFreqs)
          })



setMethod("getNucFreqs", signature(subject = "DNAStringSet"),
          function(subject) {
              nucFreqs <- colSums(letterFrequency(subject,
                                                  c("A", "C", "G", "T")))
              nucFreqs <- nucFreqs/sum(nucFreqs)
              return(nucFreqs)
          })

setMethod("getNucFreqs", signature(subject = "DNAString"),
          function(subject) {
              nucFreqs <- letterFrequency(subject, c("A", "C", "G", "T"))
              nucFreqs <- nucFreqs/sum(nucFreqs)
              return(nucFreqs)
          })

setMethod("getNucFreqs", signature(subject = "BSgenomeViews"),
          function(subject) {
              nucFreqs <- letterFrequency(subject, c("A", "C", "G", "T"))
              nucFreqs <- nucFreqs/sum(nucFreqs)
              return(nucFreqs)
          })

setMethod("getNucFreqs", signature(subject = "character"),
          function(subject) {
              if (length(subject) == 1) {
                  return(getNucFreqs(DNAString(subject)))
              } else {
                  return(getNucFreqs(DNAStringSet(subject)))
              }
          })


## validate genome input

setGeneric("validate_genome_input",
           function(genome) standardGeneric("validate_genome_input"))

setMethod("validate_genome_input", signature(genome = "FaFile"),
          function(genome) {
              return(genome)
          })

setMethod("validate_genome_input", signature(genome = "BSgenome"),
          function(genome) {
              return(genome)
          })

setMethod("validate_genome_input", signature(genome = "DNAStringSet"),
          function(genome) {
              return(genome)
          })

setMethod("validate_genome_input", signature(genome = "character"),
          function(genome) {
              return(getBSgenome(genome))
          })

setMethod("validate_genome_input", signature(genome = "character"),
          function(genome) {
              if (any(is.na(genome)))
                  stop("No genome provided")
              if (length(genome) > 1){
                  stopifnot(all(genome == genome[[1]]))
                  genome <- genome[[1]]
              }
              return(getBSgenome(genome))
          })

setMethod("validate_genome_input", signature(genome = "ANY"),
          function(genome) {
              stop("genome input must be a BSgenome, DNAStringSet, or FaFile",
                   "object or a string recognized by getBSgenome")
          })

## convert_pwm to adjust background in pwm model -------------------------------

convert_pwms <- function(pwms, bg_freqs) {
    stopifnot(inherits(pwms, "PWMatrixList"))
    lapply(pwms, convert_pwm, bg_freqs)
}

convert_pwm <- function(pwm, bg_freqs) {
    type <- pwmType(pwm)
    out <- as.matrix(pwm)
    if (type == "prob") {
        norm_mat <- matrix(bg_freqs, nrow = 4, ncol = length(pwm),
                           byrow = FALSE)
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

# make sure background is correct ----------------------------------------------
check_bg <- function(bg_freqs){
    if (length(bg_freqs) != 4)
        stop("Invalid background frequencies -- should be length 4")
    if (!all.equal(sum(bg_freqs),1) || min(bg_freqs) <= 0){
        stop("Invalid background frequencies. Should sum to 1")
    }

    if (!is.null(names(bg_freqs))){
        if (!all(names(c("A","C","G","T") %in% bg_freqs))){
            stop("Background nucleotide frequencies have names that ",
                 "don't match nucleotides! (A,C,G,T)")
        } else{
            bg_freqs <- bg_freqs[c("A","C","G","T")]
        }
    }
    return(bg_freqs)
}

get_bg <- function(bg_method, subject, genome){
    if (bg_method == "subject"){
        bg <- getNucFreqs(subject)
    } else if (bg_method == "genome"){
        if (is.null(genome))
            stop("If bg is genome, then a genome argument must ",
                 "be provided!")
        genome <- validate_genome_input(genome)
        bg <- getNucFreqs(genome)
    } else{
        bg <- rep(0.25,4)
    }
    return(bg)
}

pwm_type <- function(...){
    .Deprecated("pwmType")
    pwmType(...)
}

#' pwmType
#'
#' Determines type of PWM
#' @param pwm PWMatrix object
#' @return 'log','log2', or 'frequency' depending on type of pwm
#' @export
#' @keywords internal
#' @examples
#'
#' data(example_motifs, package = "motifmatchr")
#' pwmType(TFBSTools::toPWM(example_motifs[[1]]))
#' pwmType(TFBSTools::toPWM(example_motifs[[1]], type = "prob"))
pwmType <- function(pwm) {
  # Determine whether un-logged, natural log, or log2
  if (isTRUE(all.equal(colSums(as.matrix(pwm)), 
                       rep(1, length(pwm)),
                       check.attributes = FALSE))) {
    return("frequency")
  } else if (isTRUE(all.equal(colSums(2^(as.matrix(pwm)) *
                                      matrix(bg(pwm),
                                             byrow = FALSE,
                                             ncol = length(pwm),
                                             nrow = 4)),
                              rep(1, length(pwm)), tolerance = 10^-2,
                              check.attributes = FALSE))) {
    return("log2")
  } else if (isTRUE(all.equal(colSums(exp(as.matrix(pwm)) *
                                      matrix(bg(pwm),
                                             byrow = FALSE,
                                             ncol = length(pwm),
                                             nrow = 4)),
                              rep(1, length(pwm)), tolerance = 10^-2,
                              check.attributes = FALSE))) {
    return("log")
  } else {
    stop("Can't determine format of PWM -- should be numeric ",
         "frequency summing to 1 or log or log2 odds ratio")
  }
}

