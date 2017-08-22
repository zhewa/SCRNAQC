#!/usr/bin/env Rscript


align.rsubread <- function(demultiplex.dir, index, out.format = "BAM", out.dir = "../Alignment",
                           nthreads = 16) {
  for (i in list.files(demultiplex.dir, full.names = F)) {
    print(paste(Sys.time(), "Align sample", i))
    
    # delete results from previous run
    unlink(file.path(out.dir, i), recursive = T)
    dir.create(file.path(out.dir, i), showWarnings = FALSE, recursive = T)
    
    samples <- parse.input.files(file.path(demultiplex.dir, i))
    align(index = index,
          readfile1 = samples,
          nthreads = nthreads,
          output_format = out.format,
          output_file = file.path(out.dir, i,
                                  paste0(sub(pattern = "(.*?)\\..*$",
                                             replacement = "\\1", basename(samples)),
                                         ".", out.format)))
  }
}


report <- function(out.dir) {
  for (i in list.files(out.dir)) {
    bams <- list.files(path = file.path(out.dir, i), pattern="*.BAM$", full.names = T)
    map_prob <- propmapped(bams)
    write.table(map_prob, file.path(out.dir, i, paste(Sys.Date(), i, "alignmentstats.tab", sep="_")), sep="\t")
  }
}


build.index <- function(fa.dir, prefix) {
  out.dir <- dirname(fa.dir)
  bowtie_build(references = fa.dir, outdir = out.dir, force=T, prefix = prefix)
}


align.wrapper <- function(fastq.dir, GRCh38.index, 
                          out.format = "BAM", out.dir = "../Alignment") {
  suppressPackageStartupMessages(library(Rsubread))
  align.rsubread(fastq.dir, GRCh38.index, out.format = out.format, out.dir = out.dir)
  report(out.dir)
}

