#!/usr/bin/env Rscript

## scRNA-seq pipeline
## fastq.gz to count matrix

main <- function() {
  source("demultiplex.R")
  source("align.R")
  source("count_umi.R")
  
  # demultiplex
  bc.index.file <- "barcodes.txt"
  cut.length <- 50
  min.bc.quality <- 10
  input.dir <- "/restricted/projectnb/pulmseq/fastq/170802_NB500996_0083_AH5FC2BGX3"
  output.dir <- ".."
  umi.length <- 5
  bc.length <- 6
  fname.delimiter <- "_"
  out.folder <- "Demultiplex"
  stats.out <- "demultiplex_stats"
  out.format <- "BAM"
  align.dir <- "../Alignment"
  # alignment
  GRCh38.index <- "/restricted/projectnb/cbmhive/references/Homo_Sapiens/GRCh38/Rsubread_index/GRCh38"
  gtf.file <- "../data/gtf/Homo_sapiens.GRCh38.89.chr_ercc.gtf"
  # run pipeline
  demultiplex.wrapper(bc.index.file, input.dir, stats.out, output.dir, out.folder, min.bc.quality,
                      umi.length, bc.length, cut.length, fname.delimiter)
  fastq.dir <- file.path(output.dir, out.folder)
  align.wrapper(fastq.dir, GRCh38.index, out.format, align.dir)
  count.wrapper(alignment.dir = align.dir, gtf.file, if.bam = T, output.dir = "../Count")
}




