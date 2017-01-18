context("threshold")

data("example_motifs", package = "motifmatchr")

motif1 = motifmatchr:::convert_pwm(TFBSTools::toPWM(example_motifs[[1]]),
                                   rep(0.25,4))
motif2 = motifmatchr:::convert_pwm(TFBSTools::toPWM(example_motifs[[2]]),
                                   rep(0.25,4))


test_that("threshold matches values from TFMPvalue package",{
  m1_res1 = TFMPvalue::TFMpv2sc(motif1, 0.00005, type="PWM")
  m2_res1 = TFMPvalue::TFMpv2sc(motif2, 0.00005, type="PWM")
  m1_res2 = motifmatchr:::get_thresholds(list(motif1), rep(0.25,4), 0.00005)[1]
  m2_res2 = motifmatchr:::get_thresholds(list(motif2), rep(0.25,4), 0.00005)[1]
  expect_equal(m1_res1, m1_res2, tolerance = 0.01)
  expect_equal(m2_res1, m2_res2, tolerance = 0.01)
})


