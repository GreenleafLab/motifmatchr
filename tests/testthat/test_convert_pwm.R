context("convert_pwm")

data("example_motifs", package = "motifmatchr")

pfm <- example_motifs[[1]]
new_bg <- c(0.2, 0.3, 0.3, 0.2)
probs <- TFBSTools::Matrix(TFBSTools::toPWM(pfm, type = "prob"))

test_that("log2 PWM is rescored against new background", {
    pwm <- TFBSTools::toPWM(pfm, type = "log2")
    expect_equal(motifmatchr:::convert_pwm(pwm, new_bg),
                 log2(probs / new_bg), check.attributes = FALSE)
})

test_that("natural log PWM is rescored against new background", {
    pwm <- TFBSTools::toPWM(pfm, type = "log2")
    TFBSTools::Matrix(pwm) <- TFBSTools::Matrix(pwm) * log(2)
    expect_equal(motifmatchr:::pwmType(pwm), "log")
    expect_equal(motifmatchr:::convert_pwm(pwm, new_bg),
                 log(probs / new_bg), check.attributes = FALSE)
})

test_that("frequency PWM is converted to log2 odds vs background", {
    pwm <- TFBSTools::toPWM(pfm, type = "prob")
    expect_equal(motifmatchr:::pwmType(pwm), "frequency")
    expect_equal(motifmatchr:::convert_pwm(pwm, new_bg),
                 log2(probs / new_bg), check.attributes = FALSE)
})

test_that("PWM with non-uniform stored background is corrected", {
    pwm_bg <- c(A = 0.3, C = 0.2, G = 0.2, T = 0.3)
    pwm <- TFBSTools::toPWM(pfm, type = "log2", bg = pwm_bg)
    pwm_probs <- 2^TFBSTools::Matrix(pwm) * pwm_bg
    expect_equal(motifmatchr:::convert_pwm(pwm, new_bg),
                 log2(pwm_probs / new_bg), check.attributes = FALSE)
})

test_that("uniform background leaves log2 PWM unchanged", {
    pwm <- TFBSTools::toPWM(pfm, type = "log2")
    expect_equal(motifmatchr:::convert_pwm(pwm, rep(0.25, 4)),
                 TFBSTools::Matrix(pwm))
})

test_that("zero probabilities give an informative error, not a crash", {
    pwm <- TFBSTools::toPWM(pfm, type = "prob")
    m <- TFBSTools::Matrix(pwm)
    m[1, 1] <- 0
    m[, 1] <- m[, 1] / sum(m[, 1])
    TFBSTools::Matrix(pwm) <- m
    expect_error(motifmatchr:::convert_pwm(pwm, new_bg), "pseudocounts")
    seqs <- Biostrings::DNAStringSet("ACGTACGTACGTACGTACGT")
    expect_error(matchMotifs(pwm, seqs, bg = "even"), "pseudocounts")
})

test_that("non-finite log-odds scores give an informative error", {
    pwm <- TFBSTools::toPWM(pfm, type = "log2")
    p <- probs
    p[1, 1] <- 0
    p[, 1] <- p[, 1] / sum(p[, 1])
    TFBSTools::Matrix(pwm) <- log2(p / 0.25)
    expect_equal(motifmatchr:::pwmType(pwm), "log2")
    expect_error(motifmatchr:::convert_pwm(pwm, new_bg), "non-finite")
})
