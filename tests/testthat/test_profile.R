context("profile")
library(AlleleProfileR)

#
test_that("profile", {

  # run
  samples <- AlleleProfileR.read.folders(type = "bam")
  crispr_config <- AlleleProfileR.setup(samples = samples, genes = "files/config/example_genes.csv",
                                        index = "files/index/frag.fa", cutoff = 0.005, ignore.snp = F,
                                        cut.range = 0, ignore.single = T, cutoff.large = 25, ignore.chimeric = F,
                                        chimeric.max.insertion.size = 50, suppress.messages = T)

  AlleleProfileR.batch(crispr_config, cores=1, subset = list(c(2),c(3)))

  # check allele names
  expect_equal(AlleleProfileR.sample.distribution(crispr_config, 2, 3, plotparam = 2, plot = F)[1,"allele"], "exon.6_7insT")
  # check profile, pstop TRUE?
  expect_equal(AlleleProfileR.sample.distribution.boolean(crispr_config, 2, 3, plot = FALSE, alternate = NULL, plotparam = c("wt",'fs','snp','atg','coding','stop','pstop','sm','lg','utr','cryptic','error'))[1,'pstop'], 1)

})
