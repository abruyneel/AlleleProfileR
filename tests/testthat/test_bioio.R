context("bioio")
library(AlleleProfileR)

#
test_that("bioio", {
  # comparison with samtools faidx frag.fa chr10:40-60
  expect_equal(getreferenceC("files/index/frag.fa", "chr10", 40, 60), "CTGCCCAGCTAAGCTCCCATA")
})
