context("batch")
library(AlleleProfileR)

#
test_that("batch", {

  # run
  samples <- AlleleProfileR.read.folders(type = "bam")
  crispr_config <- AlleleProfileR.setup(samples = samples, genes = "files/config/example_genes.csv",
                                        index = "files/index/frag.fa", cutoff = 0.005, ignore.snp = F,
                                        cut.range = 0, ignore.single = T, cutoff.large = 25, ignore.chimeric = F,
                                        chimeric.max.insertion.size = 50, suppress.messages = T)

  AlleleProfileR.batch(crispr_config, cores=1, subset = list(c(1),c(1:4)))

  # no WT for any of embryo1
  expect_equal(sum(AlleleProfileR.batch.summary(crispr_config, table=T, plot=F, subset = list(c(1),c(1:4)))[,"prop"]), 0)

  # check allele names
  expect_equal(AlleleProfileR.sample.distribution(crispr_config, 1, 1, plotparam = 2, plot = F)[1,"allele"], "exon.6_14del")
  expect_equal(AlleleProfileR.sample.distribution(crispr_config, 1, 2, plotparam = 2, plot = F)[1,"allele"], "exon.2del")
  expect_equal(AlleleProfileR.sample.distribution(crispr_config, 1, 3, plotparam = 2, plot = F)[1,"allele"], "exon.20_29del")
  expect_equal(AlleleProfileR.sample.distribution(crispr_config, 1, 4, plotparam = 2, plot = F)[1,"allele"], "exon.2del.42_43insG")

})
