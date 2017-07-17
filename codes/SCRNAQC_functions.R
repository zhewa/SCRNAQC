#!/usr/bin/env Rscript

# SCRNAQC functions
# Zhe Wang
# 20170628



# ggplot publication theme
# adapted from Koundinya Desiraju
# https://rpubs.com/Koundy/71792
theme_Publication <- function(base_size=12, base_family="sans") {
  (ggthemes::theme_foundation(base_size=base_size, base_family=base_family)
   + ggplot2::theme(plot.title = element_text(face = "bold",
                                              size = rel(1), hjust = 0.5),
                    text = element_text(),
                    panel.background = element_rect(colour = NA),
                    plot.background = element_rect(colour = NA),
                    panel.border = element_rect(colour = NA),
                    axis.title = element_text(face = "bold",size = rel(1)),
                    axis.title.y = element_text(angle=90,vjust =2),
                    axis.title.x = element_text(vjust = -0.2),
                    axis.text = element_text(), 
                    axis.line = element_line(colour="black"),
                    axis.ticks = element_line(),
                    panel.grid.major = element_line(colour="#f0f0f0"),
                    panel.grid.minor = element_blank(),
                    legend.key = element_rect(colour = NA),
                    legend.position = "right",
                    legend.direction = "vertical",
                    legend.key.size= unit(0.5, "cm"),
                    legend.margin = margin(0),
                    legend.title = element_text(face="bold"),
                    plot.margin=unit(c(10,5,5,5),"mm"),
                    strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
                    strip.text = element_text(face="bold")
   ))
  
}


scale_fill_Publication <- function(...){
  ggplot2::discrete_scale("fill","Publication",manual_pal(values = cpalette), ...)
  
}


scale_colour_Publication <- function(...){
  ggplot2::discrete_scale("colour","Publication",manual_pal(values = cpalette), ...)
  
}


# read a sam file
# return a data table of the file
read.sam <- function(samfile){
  # get the maximum number of columns per row
  maxncol <- max(count.fields(samfile, sep="\t", quote="", comment.char=""))
  
  # read in SAM file
  sam <- read.table(samfile, 
                    sep="\t",
                    quote="",
                    fill=T, 
                    header=F,
                    stringsAsFactors=F,
                    na.strings=NULL,
                    comment.char="@",
                    col.names=1:maxncol)
  
  colnames(sam)[1:11] <- c("qname", "flag", "rname", "position", "mapq", "cigar",
                           "rnext", "pnext", "tlen", "seq", "qual")
  
  # convert to data.table object
  samdt <- data.table(sam[,c("qname", "rname", "position")], check.names=T)
  return (samdt)
}


# correct umi mismatch
umi.mismatch.correction <- function(samdt, current.ref, umi.max.gap, umi.edit) {
  # Add inferred_umi info to reference data table
  rdt <- samdt[rname == current.ref,]
  rdt[, c("umi", "inferred_umi") := 
        tstrsplit(qname, ":", fixed=TRUE, 
                  keep=length(tstrsplit(qname, ":", fixed=TRUE)))]
  
  # Correct UMIs with sequencing errors by looking at UMIs in surrounding region
  # These reads are considered from the same fragment if distance <= umi.edit
  # get all alignment positions
  unique.pos <- sort(unique(rdt$position))
  
  unique.pos.list <- get.adj.pos.list(unique.pos, umi.max.gap)
  
  # for each IVT fragment
  for (i in unique.pos.list) {
    all.umi.count <- sort(table(c(rdt[position %in% i, umi])))
    
    if (length(all.umi.count) > 1) {
      # temporary solution with only one iteration
      # need a recursive solution for some special cases
      sdm <- stringdistmatrix(names(all.umi.count), names(all.umi.count))
      diag(sdm) <- 100
      rownames(sdm) <- names(all.umi.count)
      colnames(sdm) <- names(all.umi.count)
      
      for (j in colnames(sdm)) {
        if (min(sdm[j,]) <= umi.edit) {
          sdm.edit.ind <- max(which(sdm[j,] <= umi.edit))
          # correct current umi j within position group i
          if (which(rownames(sdm) == j) < sdm.edit.ind) {
            rdt[umi == j & position %in% i, inferred_umi := colnames(sdm)[sdm.edit.ind]]
          }
        }
      }
    }
  }
  return (rdt)
}


get.adjacent.unique.pos <- function(unique.pos, gap) {
  adjs <- diff(unique.pos) <= gap
  ind <- which(adjs == T)
  plusone <- ind + 1
  return (unique.pos[sort(union(ind, plusone))])
}


get.adj.pos.list <- function(unique.pos, gap){
  adjs <- diff(unique.pos) <= gap
  n <- sum(adjs == F) + 1
  res <- vector("list", n) 
  
  i <- 1 # Alignment position group
  j <- 1 # index of position
  while (i <= n) {
    if (j > length(adjs)) {
      res[[i]] <- c(res[[i]], unique.pos[j])
      i <- i + 1
    } else if (adjs[j] == T) {
      res[[i]] <- c(res[[i]], unique.pos[j])
    } else if (adjs[j] == F) {
      res[[i]] <- c(res[[i]], unique.pos[j])
      i <- i + 1
    } 
    j <- j + 1
  }
  return (res)
}


get.position.with.most.reads <- function(rdt.sub) {
  res <- sort(table(rdt.sub[,position]), decreasing = T)
  return (sort(as.numeric(names(res[which(res == res[1])])))[1])
}

# correct alignment position error
alignment.position.correction <- function(rdt, pos.max.gap) {
  rdt[,inferred_pos := position]
  unique.pos <- sort(unique(rdt$position))
  
  # consider only adjacent positions with gap <= pos.max.gap
  adj.unique.pos <- get.adjacent.unique.pos(unique.pos, pos.max.gap)
  pos.group <- get.adj.pos.list(adj.unique.pos, pos.max.gap)
  # correct alignment position error
  for (i in pos.group) {
    rdt.sub <- rdt[position %in% i, ]
    unique.umi.count <- table(rdt.sub$inferred_umi)
    
    for (j in 1:length(unique.umi.count)) { 
      rdt.sub.sub <- rdt.sub[inferred_umi == names(unique.umi.count)[j], ]
      rdt[position %in% i & inferred_umi == names(unique.umi.count)[j], 
          inferred_pos := get.position.with.most.reads(rdt.sub.sub)]
    }
  }
  return (rdt)
}


get.pcr.duplicates <- function(rdt) {
  reads <- rdt[,.(inferred_umi, inferred_pos, rname)]
  unique.fragments <- unique(reads[order(reads$inferred_pos)])
  for (i in 1:nrow(unique.fragments)) {
    frag <- unique.fragments[i,]
    n <- nrow(rdt[inferred_umi==frag[[1]] & inferred_pos==frag[[2]],])
    unique.fragments[i, num := n]
  }
  return (unique.fragments)
}


num.frag.per.transcript <- function(rdt, umi.max.gap) {
  reads.dt <- unique(rdt[,.(rname, inferred_pos, inferred_umi)])
  chr <- rdt$rname[1]
  
  unique.pos <- sort(unique(reads.dt$inferred_pos))
  unique.pos.list <- get.adj.pos.list(unique.pos, umi.max.gap)
  
  num.frag <- c()
  
  for (i in seq_len(length(unique.pos.list))) {
    transcripts <- unique(reads.dt[inferred_pos %in% unique.pos.list[[i]], inferred_umi])
    for (j in transcripts) {
      num.frag <- c(num.frag, 
                    nrow(reads.dt[inferred_pos %in% unique.pos.list[[i]] & inferred_umi == j,]))
    }
  }
  res.dt <- data.table(rname = chr, transcript = paste0(chr,"_",seq_len(length(num.frag))),
                       num = num.frag)
  return (res.dt)
}


num.amplified <- function(unique.num.table) {
  return (sum(unique.num.table$num >= 2))
}


plot.num.frag.per.umi <- function(umi.table, title) {
  dt <- data.table(umi.table)
  
  breaks <- pretty(range(dt$N), n = nclass.scott(dt$N), min.n = 1)
  bwidth <- diff(breaks)[1]
  
  si <- fit.neg.binomial(umi.table)[[1]]
  m <- fit.neg.binomial(umi.table)[[2]]
  
  dt[,nbinom := dnbinom(N, size = si, mu = m)]
  
  g <- ggplot(dt, aes(N)) + 
    #geom_histogram(aes(y=..count../sum(..count..)),
    geom_histogram(aes(y=..density..),
                   binwidth=bwidth, position = "identity",
                   closed="left", boundary=0, fill="white", col="black") +
    geom_line(aes(y=nbinom), col="red") +
    #scale_x_continuous(trans="log2", name="density") +
    theme_Publication() +
    ggtitle(title) +
    xlab("Number of fragments per UMI") + ylab("Density")
  return (g)
}


fit.neg.binomial <- function(umi.table) {
  ## Try to figure out what the size (theta) and mean (mu)
  tryCatch( {
      fit = suppressWarnings(fitdistr(as.numeric(umi.table), "negative binomial", lower = 0))
      return (list(fit$estimate[1], fit$estimate[2]))
    }, error = function(f) {
      fit = suppressWarnings(fitdistr(as.numeric(umi.table), "negative binomial"))
      return (list(fit$estimate[1], fit$estimate[2]))
    }
  )
}


plot.base.fraction <- function(res.table) {

  g1 <- ggseqlogo(res.table[,umi]) + theme_Publication() + ggtitle("Original UMIs of fragments") +
    xlab("UMI position")
  g2 <- ggseqlogo(res.table[,inferred_umi]) + theme_Publication() +
    ggtitle("Inferred UMIs of fragments") + xlab("UMI position")
  return (list(g1,g2))
}


plot.num.products.per.fragment <- function(pcr.products.dt) {
  # histogram of # products per fragment
  
  #Fit negative binomial
  si <- fit.neg.binomial(pcr.products.dt$num)[[1]]
  m <- fit.neg.binomial(pcr.products.dt$num)[[2]]
  
  pcr.products.dt[,nbinom := dnbinom(num, size = si, mu = m)]
  
  pcr.products.dt[,log2.num := log2(num)]
  
  breaks <- pretty(range(pcr.products.dt$log2.num),
                   n = nclass.scott(pcr.products.dt$log2.num), min.n = 1)
  bwidth <- diff(breaks)[1]
  
  g <- ggplot(pcr.products.dt, aes(log2.num)) +
  #g <- ggplot(pcr.products.dt, aes(num)) +
    geom_histogram(aes(y=..count../sum(..count..)),
    #geom_histogram(aes(y=..density..),
                   binwidth=bwidth, closed="left", boundary=0, position="identity",
                   fill="white", col="black") +
    geom_line(aes(y=nbinom), col="red") +
    #geom_line(aes(x=log2.num, y=dnbinom(log2.num, size = si, mu = m)), col="red") +
    ggtitle("# PCR products per fragment") + 
    #scale_x_log10() +
    xlab(expression(bold(Log[2](Number~of~PCR~products)))) +
    ylab("Probability") + theme_Publication()
  
  return (g)
}


plot.num.fragments.per.transcript <- function(num.frag.table) {
  # histogram of # fragments per transcript
  
  #Fit negative binomial
  #si <- fit.neg.binomial(num.frag.table$num)[[1]]
  #m <- fit.neg.binomial(num.frag.table$num)[[2]]
  
  #num.frag.table[,nbinom := dnbinom(num, size = si, mu = m)]
  
  #num.frag.table[,log2.num := log2(num)]
  
  #breaks <- pretty(range(num.frag.table$log2.num),
  #                 n = nclass.scott(num.frag.table$log2.num), min.n = 1)
  
  breaks <- pretty(range(num.frag.table$num),
                   n = nclass.scott(num.frag.table$num), min.n = 1)
  
  bwidth <- diff(breaks)[1]
  
  #g <- ggplot(num.frag.table, aes(log2.num)) +
  g <- ggplot(num.frag.table, aes(num)) +
    geom_histogram(aes(y=..count../sum(..count..)),
    #geom_histogram(aes(y=..density..),
                   binwidth=bwidth, closed="left", boundary=0, position="identity",
                   fill="white", col="black") +
    #geom_line(aes(y=nbinom), col="red") +
    #geom_line(aes(x=log2.num, y=dnbinom(log2.num, size = si, mu = m)), col="red") +
    ggtitle("# IVT fragments per transcript") + 
    #scale_x_log10() +
    #xlab(expression(bold(Log[2](Number~of~IVT~transcripts)))) +
    xlab("Number of IVT transcripts") +
    ylab("Probability") + theme_Publication()
  
  return (g)
}


plot.stats.sam <- function(res.table, num.pcr.products.table, num.frag.table, fname) {
  # calculate fragment number
  ori.fragments <- unique(res.table[,.(umi, position)])
  inf.fragments <- unique(res.table[,.(inferred_umi, inferred_pos)])
  
  # Calculate relative abundance of UMIs
  umi.table <- table(ori.fragments$umi)
  umi.inferred.table <- table(inf.fragments$inferred_umi)
  
  nc = 2
  
  g1 <- plot.num.frag.per.umi(umi.table, "Original UMIs")
  
  g2 <- plot.num.frag.per.umi(umi.inferred.table, "Inferred UMIs")
  
  g3.g4 <- plot.base.fraction(res.table)
  
  g5 <- plot.num.products.per.fragment(num.pcr.products.table)
  
  g6 <- plot.num.fragments.per.transcript(num.frag.table)
  
  return (gridExtra::arrangeGrob(grobs = list(g1,g2,g3.g4[[1]], g3.g4[[2]], g5, g6),
                                 ncol = nc, top = grid::textGrob(paste0(fname,".sam"))))
}


get.umi.abun <- function(res.table) {
  # calculate fragment number
  ori.fragments <- unique(res.table[,.(umi, position)])
  inf.fragments <- unique(res.table[,.(inferred_umi, inferred_pos)])
  
  # Calculate relative abundance of UMIs
  umi.table <- table(ori.fragments$umi)
  umi.inferred.table <- table(inf.fragments$inferred_umi)
  return (list(umi.table, umi.inferred.table))
}


QC.sam <- function(sam, umi.edit, umi.max.gap, pos.max.gap, output.dir) {
  # read in SAM file
  samdt <- read.sam(sam)
  fname <- strsplit(last(strsplit(sam, split = "/")[[1]]), split = "\\.")[[1]][1]
  #umi.length <- as.numeric(nchar(samdt[1,tstrsplit(qname, ":", fixed=TRUE,
  #                                          keep=length(tstrsplit(qname, ":", fixed=TRUE)))]))
  # output file names
  umi.stats <- paste0(output.dir, fname, "_UMI_stats.tab")
  
  umi.qc <- paste0(output.dir, fname, "_UMI_QC_plots.pdf")
  
  
  # get reference sequence names
  chr <- mixedsort(setdiff(unique(samdt[,rname]), "*"))
  res.table <- list()
  num.pcr.products.table <- list()
  num.frag.table <- list()
  
  
  # for each reference sequence name (chr)
  for (current.ref in chr) {
    # correct umi mismatch
    rdt <- umi.mismatch.correction(samdt, current.ref, umi.max.gap, umi.edit)
    # correct alignment position error
    rdt <- alignment.position.correction(rdt, pos.max.gap)
    
    # Metrics to report:
    # Number of PCR priducts per IVT fragment
    num.pcr.products.table <- rbindlist(list(num.pcr.products.table, get.pcr.duplicates(rdt)),
                                        use.names=F, fill=F, idcol=F)
    
    # Number of IVT fragments per transcript
    num.frag.table <- rbindlist(list(num.frag.table, 
                                     num.frag.per.transcript(rdt, umi.max.gap)),
                                use.names=F, fill=F, idcol=F)
    
    # Distribution of average UMI edit distance of all reads?
    
    
    res.table <- rbindlist(list(res.table, 
                                rdt[,c("rname", "position","inferred_pos",
                                       "umi", "inferred_umi")]),
                           use.names=F, fill=F, idcol=F)
  }
  
  # Number of reads with mismatches in UMI
  num.umi.mismatch <- nrow(res.table[umi != inferred_umi,])
  
  # Number of reads with shifts in alignment position
  num.pos.shift <- nrow(res.table[position != inferred_pos,])
  
  negbinom.fit <- fit.neg.binomial(get.umi.abun(res.table)[[1]])
  negbinom.fit.inferred <- fit.neg.binomial(get.umi.abun(res.table)[[2]])
  
  prod.per.frag.fit <- fit.neg.binomial(num.pcr.products.table$num)
  
  #frag.per.trans.fit <- fit.neg.binomial(num.frag.table$num)
  
  stats.label <- c("filename",
                   "num.aligned.reads",
                   "num.unique.fragments",
                   "percent.unique.fragments",
                   "num.amplified.fragments",
                   "percent.amplified.fragments",
                   "num.unique.UMI",
                   "num.unique.UMI.corrected",
                   "size.negbinom.fit.frag.per.UMI",
                   "mu.negbinom.fit.frag.per.UMI",
                   "size.negbinom.fit.frag.per.UMI.inferred",
                   "mu.negbinom.fit.frag.per.UMI.inferred",
                   "num.reads.UMI.mismatch",
                   "percent.reads.UMI.mismatch",
                   "num.reads.pos.shift",
                   "percent.reads.pos.shift",
                   "avg.products.per.fragment",
                   "median.products.per.fragment",
                   "size.negbinom.fit.prod.per.frag.corrected",
                   "mu.negbinom.fit.prod.per.frag.corrected",
                   "num.transcripts",
                   "num.amplified.transcripts",
                   "percent.amplified.transcripts",
                   "avg.fragments.per.transcript",
                   "median.fragments.per.transcript")
                   #"size.negbinom.fit.frag.per.trans.corrected",
                   #"mu.negbinom.fit.frag.per.trans.corrected")
  
  stats <- c(fname,
             nrow(res.table),
             nrow(num.pcr.products.table),
             nrow(num.pcr.products.table)/nrow(res.table),
             num.amplified(num.pcr.products.table),
             num.amplified(num.pcr.products.table)/nrow(num.pcr.products.table),
             length(unique(res.table$umi)), 
             length(unique(res.table$inferred_umi)),
             negbinom.fit[[1]],
             negbinom.fit[[2]],
             negbinom.fit.inferred[[1]],
             negbinom.fit.inferred[[2]],
             num.umi.mismatch,
             num.umi.mismatch/nrow(res.table),
             num.pos.shift,
             num.pos.shift/nrow(res.table),
             mean(num.pcr.products.table$num),
             median(num.pcr.products.table$num),
             prod.per.frag.fit[[1]],
             prod.per.frag.fit[[2]],
             nrow(num.frag.table),
             num.amplified(num.frag.table),
             num.amplified(num.frag.table)/nrow(num.frag.table),
             mean(num.frag.table$num),
             median(num.frag.table$num))
             #frag.per.trans.fit[[1]],
             #frag.per.trans.fit[[2]]
  
  #fwrite(res.table, paste0(output.dir, last(strsplit(sam.dir, split="/")[[1]]),
  #                         "res_table.tab"), sep="\t",
  #       quote=F, row.names=F, col.names=T)
  
  dt.stats <- data.table(matrix(stats, ncol=length(stats.label), nrow=1))
  colnames(dt.stats) <- stats.label
  
  grob <- plot.stats.sam(res.table, num.pcr.products.table, num.frag.table, fname)
  return (list(dt.stats, grob))
}


# batch QC for one plate
batch.QC.sam <- function(sam.dir, umi.edit = 1, umi.max.gap = 40,
                         pos.max.gap = 5, output.dir = paste0(sam.dir,"SCRNAQC_res/")) {
  # Batch QC for one plate
  files <- mixedsort(list.files(sam.dir, full.names = T, pattern = ".sam$"))
  stats.sam <- vector(length(files), mode="list")
  plots.sam <- vector(length(files), mode="list")
  
  for (i in 1:length(files)) {
    cat("Now processing", files[i], "...\n")
    sam.res <- QC.sam(files[i], umi.edit, umi.max.gap, pos.max.gap, output.dir)
    stats.sam[[i]] <- sam.res[[1]]
    plots.sam[[i]] <- sam.res[[2]]
  }
  
  stats.res <- rbindlist(stats.sam, use.names=F, fill=F, idcol=F)
  
  fwrite(stats.res, paste0(output.dir, last(strsplit(sam.dir, split="/")[[1]]),
                           "_UMI_stats.tab"), sep="\t",
         quote=F, row.names=F, col.names=T)
  
  pdf(paste0(output.dir, last(strsplit(sam.dir, split="/")[[1]]), ".pdf"))
  for (i in 1:length(files)) {
    grid.arrange(plots.sam[[i]])
  }
  dev.off()
  
  cat("QC Done!\n")
}


visualize.QC.stats <- function(stats.file, platename) {
  dt <- fread(stats.file)
  #colnames(dt)
  
  g1 <- ggplot(dt, aes(1:nrow(dt), num.aligned.reads)) +
    geom_point() + theme_Publication() + scale_y_continuous(labels = comma) +
    ggtitle("Number of aligned reads") +
    xlab("Cell index") + ylab("Number")
  
  g2 <- ggplot(dt, aes(1:nrow(dt), percent.unique.fragments)) +
    geom_point() + geom_abline(slope=0, intercept=median(dt$percent.unique.fragments),
                               color='#E41A1C') +
    theme_Publication() + scale_y_continuous(labels = comma) +
    ggtitle("Fraction of unique fragments") +
    xlab("Cell index") + ylab("Fraction")
  
  g3 <- ggplot(dt, aes(1:nrow(dt), percent.amplified.fragments)) +
    geom_point() + geom_abline(slope=0, intercept=median(dt$percent.amplified.fragments),
                               color='#E41A1C') +
    theme_Publication() + scale_y_continuous(labels = comma) +
    ggtitle("Fraction of PCR amplified fragments") +
    xlab("Cell index") + ylab("Fraction")
  
  g4 <- ggplot(dt, aes(1:nrow(dt), num.unique.UMI.corrected)) +
    geom_point() + geom_abline(slope=0, intercept=median(dt$num.unique.UMI.corrected),
                               color='#E41A1C') +
    theme_Publication() + scale_y_continuous(labels = comma) +
    ggtitle("Number of unique UMIs after correction") +
    xlab("Cell index") + ylab("Number")
  
  g5 <- ggplot(dt, aes(1:nrow(dt), size.negbinom.fit.frag.per.UMI.inferred)) +
    geom_point() + geom_abline(slope=0,
                               intercept=median(dt$size.negbinom.fit.frag.per.UMI.inferred),
                               color='#E41A1C') +
    theme_Publication() + scale_y_continuous(labels = comma) +
    ggtitle(expression(bold(theta~of~negative~binomial~fit~to~number~of~fragments~per~UMI))) +
    xlab("Cell index") + ylab(expression(bold(theta)))
  
  g6 <- ggplot(dt, aes(1:nrow(dt), mu.negbinom.fit.frag.per.UMI.inferred)) +
    geom_point() + geom_abline(slope=0,
                               intercept=median(dt$mu.negbinom.fit.frag.per.UMI.inferred),
                               color='#E41A1C') +
    theme_Publication() + scale_y_continuous(labels = comma) +
    ggtitle(expression(bold(mu~of~negative~binomial~fit~to~number~of~fragments~per~UMI))) +
    #ggtitle("Mu of negative binomial fit to # fragments per UMI") +
    xlab("Cell index") + ylab(expression(bold(mu)))
  
  g7 <- ggplot(dt, aes(1:nrow(dt), percent.reads.UMI.mismatch)) +
    geom_point() + geom_abline(slope=0, intercept=median(dt$percent.reads.UMI.mismatch),
                               color='#E41A1C') +
    theme_Publication() + scale_y_continuous(labels = comma) +
    ggtitle("Fraction of reads with mismatch in UMI") +
    xlab("Cell index") + ylab("Fraction")
  
  g8 <- ggplot(dt, aes(1:nrow(dt), percent.reads.pos.shift)) +
    geom_point() + geom_abline(slope=0, intercept=median(dt$percent.reads.pos.shift),
                               color='#E41A1C') +
    theme_Publication() + scale_y_continuous(labels = comma) +
    ggtitle("Fraction of reads with shift in alignment position") +
    xlab("Cell index") + ylab("Fraction")
  
  g9 <- ggplot(dt, aes(1:nrow(dt), avg.products.per.fragment)) +
    geom_point() + geom_abline(slope=0, intercept=median(dt$avg.products.per.fragment),
                               color='#E41A1C') +
    theme_Publication() + scale_y_continuous(labels = comma) +
    ggtitle("Average number of PCR products per IVT fragment") +
    xlab("Cell index") + ylab("Average")
  
  g10 <- ggplot(dt, aes(1:nrow(dt), median.products.per.fragment)) +
    geom_point() + theme_Publication() + scale_y_continuous(labels = comma) +
    ggtitle("Median number of PCR products per IVT fragment") +
    xlab("Cell index") + ylab("Median")
  
  g11 <- ggplot(dt, aes(1:nrow(dt), num.transcripts)) +
    geom_point() + geom_abline(slope=0, intercept=median(dt$num.transcripts),
                               color='#E41A1C') +
    theme_Publication() + scale_y_continuous(labels = comma) +
    ggtitle("Number of transcripts") +
    xlab("Cell index") + ylab("Number")
  
  g12 <- ggplot(dt, aes(1:nrow(dt), percent.amplified.transcripts)) +
    geom_point() + geom_abline(slope=0, intercept=median(dt$percent.amplified.transcripts),
                               color='#E41A1C') +
    theme_Publication() + scale_y_continuous(labels = comma) +
    ggtitle("Fraction of IVT amplified transcripts") +
    xlab("Cell index") + ylab("Fraction")
  
  g13 <- ggplot(dt, aes(1:nrow(dt), avg.fragments.per.transcript)) +
    geom_point() + geom_abline(slope=0, intercept=median(dt$avg.fragments.per.transcript),
                               color='#E41A1C') +
    theme_Publication() + scale_y_continuous(labels = comma) +
    ggtitle("Average number of fragments per transcript") +
    xlab("Cell index") + ylab("Average")

  return (gridExtra::marrangeGrob(grobs = list(g1,g2,g3,g4,g5,g6,g7,g8,g9,g10,g11,g12,g13), 
                                  ncol = 1, nrow = 2, 
                                 top = grid::textGrob(platename)))
}





