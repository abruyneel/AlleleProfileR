context("HDR")
library(AlleleProfileR)

#
test_that("HDR", {
  # config
  samples <- AlleleProfileR.read.folders(type = "bam")
  crispr_config <- AlleleProfileR.setup(samples = samples, genes = "files/config/example_genes.csv",
                                        index = "files/index/frag.fa", cutoff = 0.005, ignore.snp = F,
                                        cut.range = 0, ignore.single = T, cutoff.large = 25, ignore.chimeric = F,
                                        chimeric.max.insertion.size = 50, suppress.messages = T)
  alternateinfo <- AlleleProfileR.alternatereference(crispr_config, alternate = "files/index/alternate.fastq", overwrite = F)

  # run
  AlleleProfileR.batch(crispr_config, cores=1, subset = list(c(3),c(2)))

  # test whether files were read properly
  expect_equal(alternateinfo[2,"allele"], "exon.22_24del")
  # test HDR reads
  expect_equal(AlleleProfileR.plot.mutationtypes(crispr_config, 3, 2, alternate = alternateinfo, title = F, plot = F)[3,'reads'], 38)
})
