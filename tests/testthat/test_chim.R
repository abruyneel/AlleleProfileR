context("chim")
library(AlleleProfileR)

#
test_that("chim", {

  # run
  samples <- AlleleProfileR.read.folders(type = "bam")
  crispr_config <- AlleleProfileR.setup(samples = samples, genes = "files/config/example_genes.csv",
                                        index = "files/index/frag.fa", cutoff = 0.005, ignore.snp = F,
                                        cut.range = 0, ignore.single = T, cutoff.large = 25, ignore.chimeric = F,
                                        chimeric.max.insertion.size = 50, suppress.messages = T)

  AlleleProfileR.batch(crispr_config, cores=1, subset = list(c(4),c(1)))

  # check allele names
  tmp <- AlleleProfileR.sample.distribution(crispr_config, 4, 1, plotparam = 2, plot = F)[,"allele"]
  expect_equal("exon.4_33del" %in% tmp, TRUE) # non-chimeric read
  expect_equal("exon.7_126del" %in% tmp, TRUE) # chimeric read

})
