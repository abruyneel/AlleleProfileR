# AlleleProfileR

Gene editing strategies, such as zinc-finger nucleases (ZFNs), transcription activator-like effector nucleases (TALENs) and clustered regularly interspaced short palindromic repeat/Cas9 (CRISPR/Cas9), are revolutionizing biochemistry. Potent, yet flexible, bioinformatics technologies are required to aid in the profiling of the generated transgenic cell lines, plants or animals. Here, we present AlleleProfileR, a novel analysis tool, written in a combination of R and C++, with the ability to batch process the sequence analysis of large and complex experiments, including base editing and large deletions created by using two guides.

## Get started
Currently, AlleleProfileR is only available from GitHub. To install it using devtools in R:
```
# install dependencies
install.packages(c("devtools","BiocManager"))
BiocManager::install()
BiocManager::install(c("BiocGenerics", "Biostrings", "GenomicAlignments", "GenomicRanges", "Rsamtools", "XVector"))
# install AlleleProfileR using devtools and github
devtools::install_github("abruyneel/AlleleProfileR")
```
This package has several dependencies, including other R-packages (such as Bioconductor: https://bioconductor.org) as well as external software. For analysing sequencing data, samtools (http://www.htslib.org), bwa (http://bio-bwa.sourceforge.net/bwa.shtml), pear (http://www.exelixis-lab.org/web/software/pear), and fastp (https://github.com/OpenGene/fastp), are needed (or similar tools). In addition, if you would like to conduct in silico experiments, wgsim (https://github.com/lh3/wgsim) is also required.

## Configuration

### File structure
AlleleProfileR requires a particular file structure for operation. In the working directory, a folder 'files' with 'config', 'index', 'input', 'output' as subfolders is required. The config folder must contain a csv file with details on the genes of interest (Table below). Reads are selected from the .bam-file for analysis if they span the region of interest as defined by the 'Start' and 'Stop' locations in the genes table. To assess whether the coding frame is altered, information on the theoretic start and stop codons of the exon need to be reported. 'StartType' or 'StopType' equals to N for normal codons: ATG, and TAG, TGA, and TAA respectively. 'PCRStart' and 'PCRStop' are the locations of the PCR primers used for amplification. 'CutSites' denotes the predicted cut positions of the editing strategy. If there are multiple cut sites, they should be separated by ';'.

| Gene  | Chr | Start | Stop | ATG | StopCodon | StartType | StartShift |StopType | StopShift | PCRStart | PCRStop | CutSites |
| ------ | ------ |------ |------ |------ |------ | ------ | ------- | ------- | ------ | ------ | ------ | ------ |
| PLN  | chr10  | 6150 | 6350  | 6179  | 6335  | N  | 0 | N  | 0 | 5500  | 6700 | NA |


The index folder must contain the reference genome in fasta format. It is not necessary to have the entire reference genome available. It is sufficient to provide the relevant chromosomes or the relevant regions within those chromosomes. The input folder must contain subfolders per sample comprising the fastq files or the bam files. The AlleleProfileR.read.folders() command scans the files/input folder for subfolders and .fastq.gz (type = "fastq") or .bam (type = "bam") files. The output folder will contain all the generated output, such as data tables and log files, but preprocessed fastq files and generated bam files will be saved into the sample input folders.

### Configuration object
The analysis functions obtain their settings from a configuration list object, set using the AlleleProfileR.setup() command. By default, the sequence analysis scope of AlleProfileR is limited to the Start and Stop positions as defined in the gene configuration. The scope can be further limited by setting the cut.range value. This will restrict AlleleProfileR's analysis to this range around the cut site(s). Morover, a cutoff value (percent value, 0 to 1) can also be set for calling variants, such that infrequent variants or sequencing noise will be ignored. Variants occuring only once can be ignored by setting the ignore.single parameter to TRUE, SNPs will be ignored by setting ignore.snp to TRUE, and chimeric reads will be ignored by setting ignore.chimeric to TRUE. Finally, a list of alternate objects can be supplied to the alternate parameter to determine proportion homology-directed repair (HDR) in gene editing experiments.

## Basic analysis example
The basic analysis operation comprises following steps: (1) pre-process the reads, (2) align reads to reference, (3) determine allelic variants, and (4) plot the results and generate summary statistics. Step (1) and (2) can be completed using external software, or by using AlleleProfileR.

### Process reads: merge and align to index
Typical read lengths are short (75-150 bp). However, when using two guides it is feasible to delete sequences larger than the read length. Additionally, the exon and genetic region of interest may be larger than the length of a typical read. If overlapping paired-end reads are merged first, the reads can be utilized for the detection of larger indels. In the current working example, we simulated paired-end fastq reads and utilized PEAR for merging the reads. Alternative algorithms may also be applied, but AlleleProfileR needs to be provided with single end or merged paired-end reads. 

AlleleProfileR utilizes reference aligned reads generated by BWA, rather than local alignments. PEAR and BWA commands can be executed on the sample datasets directly from R, either on the entire files/input folder by setting the parameter subset to NULL, or on a selection of samples by setting subset to the indices of the samples (samplestable).

```
# load library
library(AlleleProfileR)

# this does every step at once
samples <- AlleleProfileR.read.folders(type = "fastq")
AlleleProfileR.preprocess(samples, index = "files/index/frag.fa", 
                          method.qc = "fastp", params.qc = "",
                          method.merge = "pear", params.merge = "-v 30",
                          method.map = "bwa", subset = NULL)
                          
# set configuration
samplestable <- AlleleProfileR.read.folders(type = "bam")
crispr_config <- AlleleProfileR.setup(samples = samplestable, 
                                      genes = "files/config/example_genes.csv",
                                      index = "files/index/frag.fa", 
                                      cutoff = 0, 
                                      ignore.snp = F,
                                      cut.range = 0, 
                                      ignore.single = T, 
                                      cutoff.large = 25,
                                      chimeric.max.insertion.size = 250, 
                                      suppress.messages = F)
```

### Determine allelic variants
The AlleleProfileR.batch() function executes the variant determination algorithm on all data as set by the configuration or its parameters. The computation speed can be enhanced by increasing the cores value, which will split the task across multiple CPUs. The subset parameter allows the user to alter the queued files for analysis: all samples will be analysed if subset is set to NULL. Alternatively, to analyse a subset of samples for a selection of genes, subset should be set to a list containing a vector with the indices of the samples as first element, and a vector containing the indices of the genes of interest as second element. 

```
# process files and determine allelic variants
AlleleProfileR.batch(crispr_config, cores=3, subset = NULL)
```

### Summaries and plots
The variant calling algorithm writes all output to .csv-files which can be imported into R or any other software for further processing or plotting. AlleleProfileR contains several tools to visualize the results, as discussed in the sections below. 

The AlleleProfileR.batch.summary() command plots an overview of the percentage WT sequence per embryo and per gene (default: param = "wt"), and lists the number of alleles detected. Allele distributions, characterisation of the alleles and alignments for individual samples and genes can be plotted using AlleleProfileR.plot() command. Allele characterization is reported using a boolean tile plot: blue indicates TRUE, white indicates FALSE, and gray indicates NA. Alignments can be plotted using AlleleProfileR.plot.alignment(), whereas AlleleProfileR.sample.distribution() and AlleleProfileR.sample.distribution.boolean() will plot only the distributions or characterization, respectively. Finally, the read depth for a certain characterization can be plotted using AlleleProfileR.plot.readdepth().

## Reference

Bruyneel et al., "AlleleProfileR: a Versatile Tool to Identify and Profile Sequence Variants in Edited Genomes", paper submitted.
