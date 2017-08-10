#!/usr/bin/env Rscript

# count unique UMI:transcript pairs for scRNA-seq alignment files
# outputs expression matrix with columns as the cells and rows as the transcript IDs
# Zhe Wang
# 20170809


# read gtf database and return feature GRangesList by gene ID
gtf.db.read <- function(gtf.file, format="auto") {
  if (!(endsWith(gtf.file, ".gtf"))) {
    stop("Filename must end with '.gtf' or '.gff'")
  }
  gtf.db.file <- paste0(substr(gtf.file, 1, nchar(gtf.file)-3), "sqlite")
  if ((!(file.exists(gtf.file))) & (!(file.exists(gtf.db.file)))) {
    stop(paste("File", gtf.file, "does not exist"))
  }
  if (!(file.exists(gtf.db.file))) {
    print(paste("Database file", gtf.db.file, "does not exist"))
    print(paste("Generating database file", gtf.db.file))
    gtf.db <- makeTxDbFromGFF(file=gtf.file, format=format)
    saveDb(gtf.db, file=gtf.db.file)
    return (exonsBy(gtf.db, by="gene"))
  }
  gtf.db <- tryCatch(loadDb(gtf.db.file),
            error=function(e) stop(paste("Error loading database file.
Delete the file", gtf.db.file, "and try again.")))
  return (exonsBy(gtf.db, by="gene"))
}


convert.to.bam <- function(sam.list, overwrite=F, index=T) {
  for (i in sam.list) {
    tryCatch(asBam(i, overwrite=overwrite, indexDestination=index),
             error=function(e) {} )
  }
}


count.umi <- function(alignment.folder.dir, features, if.bam=F, output) {
  if (!(if.bam)) {
    sams <- mixedsort(list.files(alignment.folder.dir, full.names=T,
                                 pattern=("(.sam|.SAM)$")))
    convert.to.bam(sams)
  }
  
  bams <- mixedsort(list.files(alignment.folder.dir, full.names=T,
                               pattern=("(.bam|.BAM)$")))
  for (i in seq_len(length(bams))) {
    print(paste("Processing file", bams[i]))
    bfl <- BamFile(bams[i])
    bamGA <- readGAlignments(bfl, use.names=T)
    names(bamGA) <- data.table::last(tstrsplit(names(bamGA), ":"))
    ol = findOverlaps(features, bamGA)
    ol.dt <- data.table(gene.id=names(features)[queryHits(ol)],
                         umi=names(bamGA)[subjectHits(ol)],
                         pos=start(bamGA)[subjectHits(ol)],
                         hits=subjectHits(ol))

    # remove ambiguous gene alignments
    ol.dt <- ol.dt[!(duplicated(ol.dt, by="hits") |
                       duplicated(ol.dt, by="hits", fromLast = TRUE)), ]
    count.umi <- table(unique(ol.dt[,.(gene.id, umi)])[,gene.id])
    
    if (i == 1) {
      count.umi.dt <- data.table(gene.id=names(features))
    }
    
    count.umi.dt[[sub(pattern = "(.*?)\\..*$",
                          replacement = "\\1",
                          basename(bams[i]))]] <- 0
    count.umi.dt[gene.id %in% names(count.umi),
                     eval(sub(pattern = "(.*?)\\..*$",
                              replacement = "\\1",
                              basename(bams[i]))) := as.numeric(count.umi[gene.id])]

  }
  fwrite(count.umi.dt, output, sep="\t")
  print("Umi counting done!")
  return (count.umi.dt)
}


main <- function() {
  suppressPackageStartupMessages(library(data.table))
  suppressPackageStartupMessages(library(GenomicAlignments))
  suppressPackageStartupMessages(library(GenomicFeatures))
  suppressPackageStartupMessages(library(Rsamtools))
  suppressPackageStartupMessages(library(gtools))
  gtf.file <- "../data/gtf/Homo_sapiens.GRCh38.89.chr_ercc.gtf"
  alignment.folder.dir <- "../data/bam"
  
  features <- gtf.db.read(gtf.file)
  res <- count.umi(alignment.folder.dir, features, if.bam=F,
                   output="../res/count/0802_RPI1.tab")
}
