#!/usr/bin/env Rscript

# demultiplex fastq.gz 
# 
# Zhe Wang
# 20170811


demultiplex <- function(bc_index_file, input_dir, stats_out, output_dir, min_bc_quality,
                        umi_length=0, bc_length=6, cut_length=50, fname_delimiter="_") {
  input_files <- mixedsort(list.files(input_dir, full.names = T,
              pattern="^[^Undetermined].*.fastq$|^[^Undetermined].*.fastq.gz$", ignore.case=T))
  barcode.dt <- fread(bc_index_file, col.names = c("cell_num", "barcode"))
  meta.dt <- c()
  for (i in seq_len(length(input_files))) {
    meta.dt <- rbindlist(list(meta.dt, parse.fname(input_files[i])),
                      use.names=T, fill=F)
  }
  meta.dt <- meta.dt[order(id),]
  
  for (i in meta.dt[,unique(id)]) {
    sample.meta.dt <- meta.dt[id==i]
    lanes <- unique(sample.meta.dt[,lane])
    total <- 0
    unqualified <- 0
    
    
    for (j in lanes) {
      f1 <- sample.meta.dt[lane==j & read=="R1", fname]
      f2 <- sample.meta.dt[lane==j & read=="R2", fname]
      fq1 <- FastqStreamer(f1)
      fq2 <- FastqStreamer(f2)
      repeat {
        fqy1 <- yield(fq1)
        fqy2 <- yield(fq2)
        if ((length(fqy1) == 0 & length(fqy2) != 0) | (length(fqy1) != 0 & length(fqy2) == 0))
          stop(paste("Unequal read lengths between read1 and read2 fastq files:", f1, f2))
        else if (length(fqy1) == 0 & length(fqy2) == 0)
          break
        
        total <- total + length(fqy1)
        
        
        min.base.phred1 <- min(as(PhredQuality(substr(fqy1@quality@quality,
                                                      1, umi_length+bc_length)), "IntegerList"))
        
        fqy.dt <- data.table(rname1=tstrsplit(fqy1@id, " ")[[1]],
                             rname2=tstrsplit(fqy2@id, " ")[[1]],
                             umi=substr(fqy1@sread, 1, umi_length),
                             barcode=substr(fqy1@sread, umi_length+1, umi_length+bc_length),
                             min.phred1=min.base.phred1,
                             length1=width(fqy1),
                             read2=as.character(fqy2@sread),
                             qtring2=as.character(fqy2@quality@quality),
                             length2=width(fqy2))
        
        
      
        if (!(all(fqy.dt[, rname1] == fqy.dt[, rname2])))
          stop(paste("Abort. Read1 and read2 have different ids in files:", f1, "and", f2))
        
      
      
    }
    

    
      
      
      
      
      
    }
  }
  
  
  
  
  
  fq <- readFastq(input_files)
  fs <- FastqStreamer(input_files, verbose=T)
  fqy <- yield(fs)
}


# parse fastq filenames
# extract project name, sample ID, sample number, lane, read
# fastq files have specific naming convention
# project name, sample ID delimiter: "-" or none
# other fields delimiter: "_"
# sample fastq names:
# GD-0802-04_S4_L002_R1_001.fastq.gz
# GD-0802-04_S4_L002_R2_001.fastq.gz
# TEST3_S3_R1_001.fastq.gz
# TEST3_S3_R2_001.fastq.gz
parse.fname <- function(fastq_filename) {
  fname <- sub(pattern = "(.*?)\\..*$", replacement = "\\1", basename(fastq_filename))
  fsplit <- strsplit(fname, "_")[[1]]
  if (length(fsplit) == 5) {
    pr <- strsplit(fsplit[1], "-")[[1]]
    project <- paste(head(pr, length(pr)-1), collapse="-")
    id <- tail(pr, 1)
    num <- fsplit[2]
    lane <- fsplit[3]
    read <- fsplit[4]
  } else if (length(fsplit) == 4) {
    pr <- strsplit(gsub("([0-9]+)", "~\\1~", fsplit[1]), "~")[[1]]
    project <- paste(head(pr, length(pr)-1), collapse="-")
    id <- tail(pr, 1)
    num <- fsplit[2]
    lane <- NA
    read <- fsplit[3]
  } else {
    stop(paste("fastq filename error:", fastq_filename))
  }
  return (data.table(project=project, id=id, num=num, lane=lane, read=read, fname=fastq_filename))
}


parse.read1 <- function() {
  
}



# check minimal length
# if len(read1.qual) < umibc:
#  sample_counter['unqualified'] +=1

# if min(quals) >= int(min_bc_quality):
#  ### trim read to cut_length
#  if len(read2)>cut_length:
#  read2 = read2[0:cut_length]
# else 
# unqualified + 1


# if read1 barcode !%in% bc_index_file
# undetermined_R1.fastq + read1
# undetermined_R2.fastq + read2
# undetermined + 1


main <- function() {
  suppressPackageStartupMessages(library(data.table))
  suppressPackageStartupMessages(library(ShortRead))
  suppressPackageStartupMessages(library(gtools))
  bc_index_file <- "../data/fastq/barcodes.txt"
  fastq_filename <- "GD-0802-01_S1_L001_R1_001.fastq.gz"
  fastq_filename <- "TEST3_S3_R1_001.fastq.gz"
  input_files <- "../data/fastq/GD-0802-01_S1_L001_R1_001.fastq.gz"
  input_dir <- "../data/fastq"
  umi_length <- 5
  bc_length <- 6
}

