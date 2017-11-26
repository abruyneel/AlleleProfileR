# Copyright (C) 2017-2018  Arne A.N. Bruyneel
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.


#' @title Preprocess fastq.gz files
#' @description This function executes QC (fastp), merges paired reads (pear), and maps the reads to the reference (BWA). This external software needs to be installed manually and needs to be addressable from the terminal.
#' @param method.qc Method for conducting QC on raw reads. Currently only fastp is implemented. Use 'none' to skip this step.
#' @param params.qc Additional parameters for conducting QC on raw reads.
#' @param method.merge Method for merging paired reads. Currently only pear and fastq-join are implemented.
#' @param params.merge Additional parameters for merging paired reads
#' @param method.map Method for mapping reads to reference. Currently only BWA is implemented.
#' @inheritParams AlleleProfileR.map
#' @author Arne Bruyneel
#' @export
AlleleProfileR.preprocess <- function(files, index, method.qc = "fastp", params.qc = "",
                                      method.merge = "pear", params.merge = "-v 10",
                                      method.map = "bwa", subset = NULL) {

  # is the index generated?
  if(!file.exists(paste(index, ".fai", sep=""))) {
    # samtools faidxand bwa index
    cmd1 <- paste0("samtools faidx ", index)
    cmd2 <- paste0("bwa index ", index)
    message(cmd1, "\n");
    system(cmd1)
    message(cmd2, "\n");
    system(cmd2)
  }

  if (method.qc == "fastp") {
    AlleleProfileR.qc(files = files, method = method.qc, params = params.qc, subset = subset)
    files <- AlleleProfileR.read.folders(type = "fastq-clean")

  } else if (method.qc == "none") {
    # do nothing.

  } else {
    stop("Provide a suitable QC method")
  }

  if (method.merge %in% c("pear","fastq-join","flash")) {
    AlleleProfileR.merge(files = files, method = method.merge, params = params.merge, subset = subset)

  } else if (method.merge == "none") {
    # do nothing.

  } else {
    stop("Provide a suitable PE merging method")
  }

  if (method.map == "bwa") {
    AlleleProfileR.map(files = files, method = "bwa", index = index, subset = subset)

  } else if (method.map == "none") {
    # do nothing.

  } else {
    stop("Provide a suitable PE merging method")
  }

}

#' @title QC reads using external software
#' @description QC reads using external software. Currently only fastp is available. This external software needs to be installed manually and needs to be addressable from the terminal.
#' @param method QC method. Default is fastp.
#' @inheritParams AlleleProfileR.merge
#' @inheritParams AlleleProfileR.preprocess
#' @author Arne Bruyneel
#' @references
#' \itemize{
#'   \item Shifu Chen, Yanqing Zhou, Yaru Chen, Jia Gu; fastp: an ultra-fast all-in-one FASTQ preprocessor, Bioinformatics, Volume 34, Issue 17, 1 September 2018, Pages i884–i890, https://doi.org/10.1093/bioinformatics/bty560
#' }
#' @export
AlleleProfileR.qc <- function(files, method = "fastp", params = "", subset = NULL) {
  # subset
  if (is.null(subset)) {
    samplestable <- files
  } else {
    samplestable <- files[subset,]
  }

  if (method == "fastp") {
    # loop
    # basic command is fastp -i in.R1.fq.gz -I in.R2.fq.gz -o out.R1.fq.gz -O out.R2.fq.gz
    for(i in 1:dim(samplestable)[1]) {
      if (!is.na(samplestable$R2[i])) {
        cmd <- paste0("fastp ",params," -i ", samplestable$R1[i]," -I ", samplestable$R2[i]," -o files/input/", samplestable$Sample[i],"/", samplestable$Sample[i],".R1.clean.fastq.gz -O files/input/", samplestable$Sample[i],"/", samplestable$Sample[i],".R2.clean.fastq.gz")
      } else {
        cmd <- paste0("fastp ",params," -i ", samplestable$R1[i]," -o files/input/", samplestable$Sample[i],"/", samplestable$Sample[i],".R1.clean.fastq.gz")
      }

      message(cmd, "\n");
      system(cmd)
    }

  } else {
    stop("provide suitable qc method.")

  }

}

#' @title Merge paired-end reads
#' @description Merge peared end reads using PEAR, FLASH or fastq-join. This external software needs to be installed manually and needs to be addressable from the terminal.
#' @param params String with additional parameters
#' @param method Merge method, options are pear, fastq-join or flash. Default is pear (PEAR).
#' @param subset Vector with the indices of the entries that should be processed. Default is NULL, all entries will be analyzed.
#' @inheritParams AlleleProfileR.map
#' @author Arne Bruyneel
#' @references
#' \itemize{
#'   \item Zhang, J., et al. PEAR: a fast and accurate Illumina Paired-End reAd mergeR. Bioinformatics 2014;30(5):614-620.
#'   \item Magoč, T. and Salzberg, S.L. FLASH: fast length adjustment of short reads to improve genome assemblies. Bioinformatics 2011;27(21):2957-2963.
#'   \item Aronesty, E. Comparison of Sequencing Utility Programs. The Open Bioinformatics Journal  2013;7:1-8.
#' }
#' @export
AlleleProfileR.merge <- function(files, method = "pear", params = "-v 10", subset = NULL) {
  # subset
  if (is.null(subset)) {
    samplestable <- files
  } else {
    samplestable <- files[subset,]
  }

  if (method %in% c("pear","fastq-join","flash")) {
    # loop
    for(i in 1:dim(samplestable)[1]) {
      if (!is.na(samplestable$R2[i]) & !is.na(samplestable$R1[i])) {
        if (method == "pear") {
          # The basic command to run PEAR is
          # ./pear -f forward_read.fastq -r reverse_read.fastq -o output_prefix
          # the output file of interest is: output_prefix.assembled.fastq
          cmd <- paste0("pear ",params," -f ", samplestable$R1[i]," -r ", samplestable$R2[i]," -o files/input/", samplestable$Sample[i], "/", samplestable$Sample[i])
          # execute
          message(cmd, "\n");
          system(cmd)

        } else if (method == "fastq-join") {
          cmd1 <- paste0("fastq-join ",params," ", samplestable$R1[i]," ", samplestable$R2[i]," -o files/input/", samplestable$Sample[i], "/", samplestable$Sample[i], ".assembled.fastq")
          cmd2 <- paste0("cp"," files/input/", samplestable$Sample[i], "/", samplestable$Sample[i], ".assembled.fastqjoin"," files/input/", samplestable$Sample[i], "/", samplestable$Sample[i], ".assembled.fastq")
          # execute
          message(cmd1, "\n");
          system(cmd1)
          message(cmd2, "\n");
          system(cmd2)

        } else if (method == "flash") {
          cmd1 <- paste0("flash ",params, " -d files/input/", samplestable$Sample[i], "/ ", samplestable$R1[i]," ", samplestable$R2[i])
          cmd2 <- paste0("cp"," files/input/", samplestable$Sample[i], "/out.extendedFrags.fastq"," files/input/", samplestable$Sample[i], "/", samplestable$Sample[i], ".assembled.fastq")
          # execute
          message(cmd1, "\n");
          system(cmd1)
          message(cmd2, "\n");
          system(cmd2)

        }

      }
    }

  } else {
    stop("provide suitable merge method.")

  }

}

#' @title Align reads to the reference genome
#' @description Align the reads to the reference using BWA. This external software needs to be installed manually and needs to be addressable from the terminal.
#' @param files Vector with files. Output of AlleleProfileR.read.folders
#' @param method Mapping method, default is BWA (BWA-MEM).
#' @inheritParams AlleleProfileR.setup
#' @inheritParams AlleleProfileR.merge
#' @author Arne Bruyneel
#' @references
#' \itemize{
#'   \item Li, H. and Durbin, R. Fast and accurate long-read alignment with Burrows-Wheeler transform. Bioinformatics 2010;26(5):589-595.
#' }
#' @export
AlleleProfileR.map <- function(files, method = "bwa", index, subset = NULL) {
  # subset
  if (is.null(subset)) {
    samplestable <- files
  } else {
    samplestable <- files[subset,]
  }

  # index
  bwa_index <- index

  if (method == "bwa") {
    # loop
    for(i in 1:dim(samplestable)[1]) {
      # the basic command is: bwa mem ref.fa reads.fq > aln-se.sam
      cmd <- paste0("bwa mem ", bwa_index, " files/input/", samplestable$Sample[i], "/", samplestable$Sample[i], ".assembled.fastq | samtools view -Sb - > files/input/", samplestable$Sample[i], "/", samplestable$Sample[i],".bam")
      # execute
      message(cmd, "\n");
      system(cmd)
    }

  } else {
    stop("provide suitable mapping method.")

  }

}

#' @title Simulate paired-end reads from fasta files
#' @description This function samples reads from one or more fasta file using wgsim. This external software needs to be installed manually and needs to be addressable from the terminal.
#' No data is returned: fastq.gz files are generated and deposited in the input.folder directory.
#' @param input.folder Path to input fasta file folder
#' @param input.file Vector with one or more .fas file names
#' @param params Parameters for wgsim, default is: -e 0 -d 250 -N 50000 -1 150 -2 150 -r 0 -R 0 -X 0
#' @references
#' \itemize{
#'   \item https://github.com/lh3/wgsim
#' }
#' @export
AlleleProfileR.simulate <- function(input.folder, input.file, params = "-e 0 -d 250 -N 50000 -1 150 -2 150 -r 0 -R 0 -X 0") {
  # delete the data from a previous simulation
  system(paste0("rm ",input.folder,"/R1.fastq.gz",sep=""))
  system(paste0("rm ",input.folder,"/R2.fastq.gz",sep=""))

  for (i in 1:length(input.file)) {
    # deletion the existing fastq.gz files
    system(paste0("rm ",input.folder,"/",i,"_R1.fastq.gz",sep=""))
    system(paste0("rm ",input.folder,"/",i,"_R2.fastq.gz",sep=""))

    # the basic command is: bwa mem ref.fa reads.fq > aln-se.sam
    cmd <- paste0("wgsim ", params, " ", input.folder, "/", input.file[i] ," ", input.folder, "/",i,"_R1.fastq ", input.folder, "/",i,"_R2.fastq",sep="")
    message(cmd, "\n");
    system(cmd)
  }

  # for merging files:
  system(paste0("cat ", input.folder, "/*_R1.fastq > ", input.folder, "/R1.fastq", sep=""))
  system(paste0("cat ", input.folder, "/*_R2.fastq > ", input.folder, "/R2.fastq", sep=""))

  # compress
  system(paste0("gzip ", input.folder, "/R1.fastq ", input.folder, "/R1.fastq.gz", sep=""))
  system(paste0("gzip ", input.folder, "/R2.fastq ", input.folder, "/R2.fastq.gz", sep=""))

}

