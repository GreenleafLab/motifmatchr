context("match_motifs")

data("example_motifs", package = "motifmatchr")

peaks <- GenomicRanges::GRanges(seqnames = c("chr1","chr2","chr2"),
                                ranges = IRanges::IRanges(start = c(76585873,
                                                                    42772928,
                                                                    100183786),
                                                          width = 500))

motif1 = motifmatchr:::convert_pwm(TFBSTools::toPWM(example_motifs[[1]]),
                                   rep(0.25,4))
motif2 = motifmatchr:::convert_pwm(TFBSTools::toPWM(example_motifs[[2]]),
                                   rep(0.25,4))
motif3 = motifmatchr:::convert_pwm(TFBSTools::toPWM(example_motifs[[3]]),
                                   rep(0.25,4))


se <- SummarizedExperiment::SummarizedExperiment(assays =
                                                     list(counts =
                                                              matrix(1,
                                                                     ncol = 4,
                                                                     nrow = 3)),
                                                 rowRanges = peaks)

dss <-
    Biostrings::getSeq(BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19,
                          peaks)

ch <- as.character(dss)

thresh <-  motifmatchr:::get_thresholds(list(motif1, motif2, motif3),
                                        rep(0.25,4), 0.00005)[1:3]

example_pwms <- do.call(TFBSTools::PWMatrixList,lapply(example_motifs,
                                                       TFBSTools::toPWM))

bs_method <- function(motif, s, score){
  forward_matches <- Biostrings::matchPWM(motif, s, min.score = score)
  reverse_matches <- Biostrings::matchPWM(motif,
                                          Biostrings::reverseComplement(s),
                                          min.score = score)
  (length(forward_matches) !=0) || (length(reverse_matches) != 0)
}

m1 <- sapply(dss, function(x) bs_method(motif1,x, thresh[1]))
m2 <- sapply(dss, function(x) bs_method(motif2,x, thresh[2]))
m3 <- sapply(dss, function(x) bs_method(motif3,x, thresh[3]))

bs_res <- cbind(m1,m2,m3)
colnames(bs_res) <- names(example_motifs)

bs_method_scores <- function(motif, s, score){
    forward_matches <- Biostrings::matchPWM(motif, s, min.score = score,
                                            with.score = TRUE)
    reverse_matches <- Biostrings::matchPWM(motif,
                                            Biostrings::reverseComplement(s),
                                            min.score = score,
                                            with.score = TRUE)
    if ((length(forward_matches) !=0) || (length(reverse_matches) != 0))
        max(c(S4Vectors::mcols(forward_matches)$score,
              S4Vectors::mcols(reverse_matches)$score)) else 0
}

m1s <- sapply(dss, function(x) bs_method_scores(motif1,x, thresh[1]))
m2s <- sapply(dss, function(x) bs_method_scores(motif2,x, thresh[2]))
m3s <- sapply(dss, function(x) bs_method_scores(motif3,x, thresh[3]))

bs_res_scores <- cbind(m1s,m2s,m3s)
colnames(bs_res_scores) <- names(example_motifs)

bs_method_counts <- function(motif, s, score){
    forward_matches <- Biostrings::matchPWM(motif, s, min.score = score,
                                            with.score = TRUE)
    reverse_matches <- Biostrings::matchPWM(motif,
                                            Biostrings::reverseComplement(s),
                                            min.score = score,
                                            with.score = TRUE)
    length(forward_matches) + length(reverse_matches)
}

m1c <- sapply(dss, function(x) bs_method_counts(motif1,x, thresh[1]))
m2c <- sapply(dss, function(x) bs_method_counts(motif2,x, thresh[2]))
m3c <- sapply(dss, function(x) bs_method_counts(motif3,x, thresh[3]))

bs_res_counts <- cbind(m1c,m2c,m3c)
colnames(bs_res_counts) <- names(example_motifs)


bs_method_positions <- function(motif, s, score, ranges){
    forward_matches <- Biostrings::matchPWM(motif, s, min.score = score,
                                            with.score = TRUE)
    reverse_matches <- Biostrings::matchPWM(motif,
                                            Biostrings::reverseComplement(s),
                                            min.score = score,
                                            with.score = TRUE)
    if ((length(forward_matches) !=0) || (length(reverse_matches) != 0)){
        franges <- as(forward_matches,"IRanges")
        rranges <- as(reverse_matches,"IRanges")
        out <- GenomicRanges::GRanges(seqnames = seqnames(ranges),
                                      ranges = IRanges::IRanges(start = c(start(ranges) + start(franges) - 1,
                                                                          end(ranges) - end(rranges) + 1),
                                                               end = c(start(ranges) + end(franges) - 1,
                                                                       end(ranges) - start(rranges) + 1)),
                                      strand = c(rep("+",length(franges)),rep("-",length(rranges))),
                                      scores = c(S4Vectors::mcols(forward_matches)$score,
                                                  S4Vectors::mcols(reverse_matches)$score))
        out
    } else{
        GenomicRanges::GRanges()
    }
}

m1p <- do.call(c,lapply(seq_along(dss), function(x) bs_method_positions(motif1,dss[[x]], thresh[1], peaks[x])))
m2p <- do.call(c,lapply(seq_along(dss), function(x) bs_method_positions(motif2,dss[[x]], thresh[2], peaks[x])))
m3p <- do.call(c,lapply(seq_along(dss), function(x) bs_method_positions(motif3,dss[[x]], thresh[3], peaks[x])))

bs_res_positions <- list(m1p,m2p,m3p)
names(bs_res_positions) <- names(example_motifs)



# Output of matches ------------------------------------------------------------

test_that("Can run match_pwm with PFMatrixList and peaks",{
  mm_res <- match_motifs(example_motifs, peaks, bg = rep(0.25,4))
  expect_equal(as.matrix(assays(mm_res)$motif_matches), bs_res)
  expect_is(mm_res, "SummarizedExperiment")
})

test_that("Can run match_pwm with PFMatrixList and SummarizedExperiment",{
  mm_res <- match_motifs(example_motifs, se, bg = rep(0.25,4))
  expect_equal(as.matrix(assays(mm_res)$motif_matches), bs_res)
  expect_is(mm_res, "SummarizedExperiment")
})

test_that("Can run match_pwm with PFMatrixList and DNAStringSet",{
  mm_res <- match_motifs(example_motifs, dss, bg = rep(0.25,4))
  expect_equal(as.matrix(assays(mm_res)$motif_matches), bs_res)
  expect_is(mm_res, "SummarizedExperiment")
})

test_that("Can run match_pwm with PFMatrixList and character string",{
  mm_res <- match_motifs(example_motifs, ch, bg = rep(0.25,4))
  expect_equal(as.matrix(assays(mm_res)$motif_matches), bs_res)
  expect_is(mm_res, "SummarizedExperiment")
})

test_that("Can run match_pwm with PFMatrixList and DNAString",{
  mm_res <- match_motifs(example_motifs, dss[[3]], bg = rep(0.25,4))
  expect_equal(as.matrix(assays(mm_res)$motif_matches)[1,], bs_res[3,])
  expect_is(mm_res, "SummarizedExperiment")
})


test_that("Can run match_pwm with PWMatrixList and peaks",{
  mm_res <- match_motifs(example_pwms, peaks, bg = rep(0.25,4))
  expect_equal(as.matrix(assays(mm_res)$motif_matches), bs_res)
  expect_is(mm_res, "SummarizedExperiment")
})

test_that("Can run match_pwm with PWMatrixList and SummarizedExperiment",{
  mm_res <- match_motifs(example_pwms, se, bg = rep(0.25,4))
  expect_equal(as.matrix(assays(mm_res)$motif_matches), bs_res)
  expect_is(mm_res, "SummarizedExperiment")
})

test_that("Can run match_pwm with PWMatrixList and DNAStringSet",{
  mm_res <- match_motifs(example_pwms, dss, bg = rep(0.25,4))
  expect_equal(as.matrix(assays(mm_res)$motif_matches), bs_res)
  expect_is(mm_res, "SummarizedExperiment")
})

test_that("Can run match_pwm with PWMatrixList and character string",{
  mm_res <- match_motifs(example_pwms, ch, bg = rep(0.25,4))
  expect_equal(as.matrix(assays(mm_res)$motif_matches), bs_res)
})

test_that("Can run match_pwm with PWMatrixList and DNAString",{
  mm_res <- match_motifs(example_pwms, dss[[3]], bg = rep(0.25,4))
  expect_equal(as.matrix(assays(mm_res)$motif_matches)[1,], bs_res[3,])
})



test_that("Can run match_pwm with PWMatrix and peaks",{
  mm_res <- match_motifs(example_pwms[[3]], peaks, bg = rep(0.25,4))
  expect_equal(as.vector(assays(mm_res)$motif_matches), bs_res[,3])
  expect_is(mm_res, "SummarizedExperiment")
})

test_that("Can run match_pwm with PWMatrix and SummarizedExperiment",{
  mm_res <- match_motifs(example_pwms[[3]], se, bg = rep(0.25,4))
  expect_equal(as.vector(assays(mm_res)$motif_matches), bs_res[,3])
  expect_is(mm_res, "SummarizedExperiment")
})

test_that("Can run match_pwm with PWMatrix and DNAStringSet",{
  mm_res <- match_motifs(example_pwms[[3]], dss, bg = rep(0.25,4))
  expect_equal(as.vector(assays(mm_res)$motif_matches), bs_res[,3])
  expect_is(mm_res, "SummarizedExperiment")
})

test_that("Can run match_pwm with PWMatrix and character string",{
  mm_res <- match_motifs(example_pwms[[3]], ch, bg = rep(0.25,4))
  expect_equal(as.vector(assays(mm_res)$motif_matches), bs_res[,3])
  expect_is(mm_res, "SummarizedExperiment")
})

test_that("Can run match_pwm with PWMatrix and DNAString",{
  mm_res <- match_motifs(example_pwms[[3]], dss[[3]], bg = rep(0.25,4))
  expect_equal(as.vector(assays(mm_res)$motif_matches), unname(bs_res[3,3]))
})


test_that("Can run match_pwm with PFMatrix and peaks",{
  mm_res <- match_motifs(example_motifs[[3]], peaks, bg = rep(0.25,4))
  expect_equal(as.vector(assays(mm_res)$motif_matches), bs_res[,3])
  expect_is(mm_res, "SummarizedExperiment")
})

test_that("Can run match_pwm with PFMatrix and SummarizedExperiment",{
  mm_res <- match_motifs(example_motifs[[3]], se, bg = rep(0.25,4))
  expect_equal(as.vector(assays(mm_res)$motif_matches), bs_res[,3])
  expect_is(mm_res, "SummarizedExperiment")
})

test_that("Can run match_pwm with PFMatrix and DNAStringSet",{
  mm_res <- match_motifs(example_motifs[[3]], dss, bg = rep(0.25,4))
  expect_equal(as.vector(assays(mm_res)$motif_matches), bs_res[,3])
  expect_is(mm_res, "SummarizedExperiment")
})

test_that("Can run match_pwm with PFMatrix and character string",{
  mm_res <- match_motifs(example_motifs[[3]], ch, bg = rep(0.25,4))
  expect_equal(as.vector(assays(mm_res)$motif_matches), bs_res[,3])
  expect_is(mm_res, "SummarizedExperiment")
})

test_that("Can run match_pwm with PFMatrix and DNAString",{
  mm_res <- match_motifs(example_motifs[[3]], dss[[3]], bg = rep(0.25,4))
  expect_equal(as.vector(assays(mm_res)$motif_matches), unname(bs_res[3,3]))
})

# Output of scores ------------------------------------------------------------

test_that("Can run match_pwm with ouptu of scores",{
    mm_res <- match_motifs(example_motifs, peaks, bg = rep(0.25,4),
                           out = "scores")
    expect_equal(as.matrix(assays(mm_res)$motif_scores), bs_res_scores)
    expect_is(mm_res, "SummarizedExperiment")
})

# Output of positions ----------------------------------------------------------

test_that("Can run match_pwm with output of positions",{
    mm_res <- match_motifs(example_motifs, peaks, bg = rep(0.25,4),
                           out = "positions")
    expect_equal(do.call(c,lapply(mm_res, start)),
                 do.call(c,lapply(bs_res_positions, start)))
    expect_equal(do.call(c,lapply(mm_res, end)),
                 do.call(c,lapply(bs_res_positions, end)))
    expect_equal(do.call(c,lapply(mm_res, strand)),
                 do.call(c,lapply(bs_res_positions, strand)))
    expect_equal(do.call(c,lapply(mm_res, function(x) x$score)),
                 do.call(c,lapply(bs_res_positions, function(x) x$scores)))
    expect_is(mm_res, "GRangesList")
})

# Accessors --------------------------------------------------------------------

test_that("Accessors work",{
    test_obj <- SummarizedExperiment::SummarizedExperiment(assays = list(motif_matches = bs_res,
                                                                         motif_scores = bs_res_scores,
                                                                         motif_counts = bs_res_counts))
    expect_equal(motif_matches(test_obj), bs_res)
    expect_equal(motif_scores(test_obj), bs_res_scores)
    expect_equal(motif_counts(test_obj), bs_res_counts)
})


