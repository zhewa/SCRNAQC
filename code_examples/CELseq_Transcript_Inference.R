
dumi = function(files, EM="global", mer.filter=0.01, aer.filter=0.01, verbose=TRUE) {

  bfl = BamFileList(files)
  
  ## Get summary stats for all files
  bam.align = countBam(bfl, param=ScanBamParam(scanBamFlag(isUnmappedQuery=FALSE, isNotPassingQualityControls=FALSE)))
  bam.unalign = countBam(bfl, param=ScanBamParam(scanBamFlag(isUnmappedQuery=TRUE)))
  bam.low.q = countBam(bfl, param=ScanBamParam(scanBamFlag(isUnmappedQuery=FALSE, isNotPassingQualityControls=TRUE)))
  
  ## Perform EM for each cell individually

  cell.mer = rep(NA, length(files))
  cell.mer.counts = rep(NA, length(files))
  cell.aer = rep(NA, length(files))  
  cell.aer.counts = rep(NA, length(files))
    
  for(i in 1:length(files)) {
    bem = barcodeEM(files(i), verbose=verbose)
    cell.aer[i] = bem$aer
    cell.mer[i] = bem$mer
    
    
  }
}





barcodeEM = function(filename, window=10, umi.flag=5, mismatch.error.rate = 0.05, alignment.error.rate = rep(1, (window*2)+1) / ((window*2)+1), max.iter=10, verbose=TRUE) {
  suppressPackageStartupMessages(require(stringr))
  suppressPackageStartupMessages(require(GenomicAlignments))
  suppressPackageStartupMessages(require(data.table))
  suppressPackageStartupMessages(require(stringdist))

  bamGA = readGAlignments(filename, param=ScanBamParam(what=c("qname", "flag")))
  
  umi = read.umi(mcols(bamGA)$qname, umi.flag)

  position = start(bamGA)
  strand = as.factor(strand(bamGA)) == "+"
  position[!strand] = end(bamGA)[!strand]

  reads = data.table(data.frame(qname=mcols(bamGA)$qname, chr=seqnames(bamGA), strand=strand, position=position, umi=umi, inferred_position=position, inferred_umi=umi))
  reads[, `:=` (INITIAL_COUNT = .N, INITIAL_GRP=.GRP, INITIAL_GRP_IX=1:.N), by=list(chr, strand, inferred_position, inferred_umi)]
  reads[, ID := 1:.N]
           
  ## Set up shorthand variables for use in rest of script
  aer = alignment.error.rate
  names(aer) = as.character(-(window):(window))
  umi.len = nchar(umi[1])
  mer = mismatch.error.rate
  w = window
  w.size = (2*w)+1
  mer.full = list()
  aer.full = list()
  
  ## Set up variables to save info after each iteration
  previous.umi = reads$inferred_umi
  previous.position = reads$inferred_position
  
  iter = 1
  continue = TRUE
  while(continue == TRUE & iter <= max.iter) {
    
    ## Set up mismatch and alignment error probability matrix
    probs = get.probability.matrix(umi.len, mer, aer)

	## Calculate read fragment memebership based on lasted UMI/position assignment
	reads = reads[, `:=` (COUNT = .N, GRP=.GRP, GRP_IX=1:.N), by=list(chr, strand, inferred_position, inferred_umi)]
	reads.dedup = subset(reads, GRP_IX == 1)

	next.umi = rep(NA, nrow(reads.dedup))
	next.position = rep(NA, nrow(reads.dedup))
	
	## Identify most likely assignment for each read fragment group
	if(verbose) print(sprintf("Expectation ... Calculating most likely fragment memberships"))	
	for(ix in 1:nrow(reads.dedup)) {
	  current.umi = as.character(reads.dedup[ix,inferred_umi])
	  current.position = reads.dedup[ix,inferred_position]
	  current.strand = reads.dedup[ix,strand]
	  current.chr = reads.dedup[ix,chr]
	
	  reads.ix = reads.dedup$chr == current.chr & reads.dedup$inferred_position > (current.position-w) & reads.dedup$inferred_position < (current.position+w) & reads.dedup$strand == current.strand

      ## Calculate UMI edit distance and genomic distance for each read within window
      reads.dist = stringdist(current.umi, reads.dedup[reads.ix,inferred_umi])
      if(current.strand == TRUE) {
	    reads.pos = reads.dedup$inferred_position[reads.ix] - current.position
      } else {
        reads.pos = -(reads.dedup$inferred_position[reads.ix] - current.position)
      }

      reads.prob = probs[cbind(reads.dist+1, reads.pos+w+1)]
      ll = log(reads.prob*reads.dedup$COUNT[reads.ix])

	  ties = ifelse(current.strand == TRUE, "left", "right")
	  ll.max.ix = get.max.ll(ll, ties=ties)

	  next.umi[ix] = as.character(reads.dedup$inferred_umi[reads.ix][ll.max.ix])
	  next.position[ix] = reads.dedup$inferred_position[reads.ix][ll.max.ix]
	}
 
	if(verbose) print(sprintf("Maximization ... Calculating new mismatch error rate"))
	mer.temp = calcMismatchErrorRate(reads.dedup$umi, next.umi, reads.dedup$COUNT)
	mer = mer.temp[1]
    mer.full = c(mer.full, list(mer.temp))  
    
	if(verbose) print(sprintf("Maximization ... Calculating new alignment error rate"))
	aer.temp = calcAlignmentErrorRate(reads.dedup$position, next.position, reads.dedup$strand, reads.dedup$COUNT, window=w)
	aer = aer.temp[[1]]  
    aer.full = c(aer.full, list(aer.temp))
    
    id = as.numeric(reads.dedup$ID)
    reads$inferred_umi[id] = next.umi
    reads$inferred_position[id] = next.position
  
	if(verbose) print(sprintf("Completed iteration %d", iter))
	iter = iter + 1
	
	total.diff = sum(!(reads$inferred_umi == previous.umi & reads$inferred_position == previous.position))
	if(verbose) print(sprintf("Total number of reads that changed fragment membership: %d", total.diff))
	if(total.diff == 0) {
	  continue = FALSE
	}
	previous.umi = reads$inferred_umi
	previous.position = reads$inferred_position
  }
  
  ## Make a final read matrix with duplicate infomation
  diff = !(reads$umi == reads$inferred_umi & reads$position == reads$position)
  reads[,`:=` (COUNT = .N, GRP=.GRP, GRP_IX=1:.N, Is_Different=diff, Duplicate = FALSE)]
  reads = reads[,DUPLICATE_IX := 1:.N, by=list(chr, strand, inferred_position, inferred_umi, Is_Different)]
  ind = reads$DUPLICATE_IX > 1 | diff
  reads$Duplicate[ind] = TRUE
  reads[,DUPLICATE_IX := NULL]

  return(list(ga=bamGA[reads$Duplicate == FALSE], reads=reads, aer=aer, mer=mer, aer.full=aer.full, mer.full=mer.full))
}

calcMismatchErrorRate = function(actual.umi, inferred.umi, counts) {
  actual.umi = as.character(actual.umi)
  inferred.umi = as.character(inferred.umi)
  
  edit.dist = rep(0, length(actual.umi))
  ind = actual.umi != inferred.umi
  
  c = cbind(as.character(actual.umi[ind]), as.character(inferred.umi[ind]))
  edit.dist[ind] = unlist(lapply(1:nrow(c), function(i) stringDist(c(c[i,1], c[i,2]), method="hamming")[1]))
  
  mer.mismatch = sum(edit.dist*counts)
  mer.total = sum(counts*nchar(actual.umi[1]))
  new.mer =  mer.mismatch / mer.total
  
  return(c(new.mer, mer.mismatch, mer.total))
}


calcAlignmentErrorRate = function(actual.position, inferred.position, strand, counts, window=10, pseudo=1) {
  
  window.len = (2*window)+1
  new.aer = rep(0, window.len)
  names(new.aer) = (-window):(window)
  
  new.aer.forward = factor(inferred.position[strand] - actual.position[strand], levels=names(new.aer))
  new.aer.reverse = factor(-(inferred.position[!strand] - actual.position[!strand]), levels=names(new.aer))
  
  new.aer.forward.agg = aggregate(counts[strand], by=list(new.aer.forward), sum)
  new.aer.reverse.agg = aggregate(counts[!strand], by=list(new.aer.reverse), sum) 
  
  new.aer[new.aer.forward.agg[,1]] = new.aer[new.aer.forward.agg[,1]] + new.aer.forward.agg[,2]
  new.aer[new.aer.reverse.agg[,1]] = new.aer[new.aer.reverse.agg[,1]] + new.aer.reverse.agg[,2]
  
  new.aer.raw = new.aer
  new.aer = (new.aer + pseudo) / (sum(new.aer) + pseudo*window.len)
  
  return(list(new.aer, new.aer.raw))
}





get.max.ll = function(v, ties=c("left", "right")) {
  ties = match.arg(ties)

  v.max = max(v, na.rm=TRUE)
  v.max.ix = which(v == v.max)

  l = length(v.max.ix)+1
  ix = ifelse(ties == "left", floor(l/2), ceiling(l/2))
  return(v.max.ix[ix])
}



read.umi = function(qname, umi.flat) {
  if(length(umi.flag) == 1 & class(umi.flag) %in% c("numeric", "integer")) {
    umi = str_sub(qname, start=-umi.flag)
  } else if(length(umi.flag) == 2 & class(umi.flag) %in% c("numeric", "integer")) {  
    umi = str_sub(qname, start=umi.flag[1], end=umi.flag[2])
  } else if(length(umi.flag) == 1 & class(umi.flag) == "character") {
    umi = str_match(qname, pattern=umi.flag)[,2]
    if(sum(is.na(umi)) > 0) {
      stop(paste("UMIs not found for ", sum(is.na(umi)), " reads. Please check matching pattern", sep=""))
    }
  } else {
    stop("umi.flag not in a recognized format. See ?barcodeEM for details.")
  }
  return(umi)
}


get.probability.matrix = function(umi.len, mer, aer) {

  probs = do.call(cbind, lapply(1:length(aer), function(i) aer[i] * dbinom(x=0:umi.len, size=umi.len, prob=mer)))
  rownames(probs) = 0:umi.len
  colnames(probs) = names(aer)
  
  return(probs)
}

