# Obsolete stuff

umi.mismatch.correction.obsolete <- function(samdt, current.ref, umi.window, umi.edit) {
  cat("Now processing", current.ref, "...\n")
  
  # Add UMI info to reference data table
  rdt <- samdt[rname == current.ref,]
  rdt[, c("umi", "inferred_umi") := 
        tstrsplit(qname, ":", fixed=TRUE, 
                  keep=length(tstrsplit(qname, ":", fixed=TRUE)))]
  
  # Correct UMIs with sequencing errors by looking at UMIs in surrounding region
  
  # get all alignment positions
  unique.pos <- sort(unique(rdt$position))
  
  for (i in unique.pos) {
    # Get range data table for position and surrounding window
    # get all alignment positions ranging between i-umi.window to i+umi.window
    rdt.sub <- rdt[position == i, ]
    rdt.flank.5p <- rdt[position >= (i-umi.window) & position <= (i-1), ]
    rdt.flank.3p <- rdt[position >= (i+1) & position <= (i+umi.window), ]
    rdt.flank <- rbindlist(list(rdt.flank.5p, rdt.flank.3p), use.names=F, fill=F, idcol=F)
    
    # Get all unique UMIs in the region
    all.umi.count <- sort(table(c(rdt.sub[,umi], rdt.flank[,umi])))
    
    
    # Align UMIs to all other UMIs. Assign UMIs with lower counts to 
    # matching UMI with higher counts
    
    if (length(all.umi.count) > 1) {
      sdm <- stringdistmatrix(names(all.umi.count), names(all.umi.count))
      diag(sdm) <- 100
      rownames(sdm) <- names(all.umi.count)
      colnames(sdm) <- names(all.umi.count)
      
      position.umi.count <- sort(table(rdt.sub$umi))
      
      # for umi k at position i
      for (k in names(position.umi.count)) {
        
        # Get the min edit distance for that UMI
        sdm.min.edit <- min(sdm[k,])
        
        if (sdm.min.edit <= umi.edit) {
          # Get the index of the edit distance less than umi.edit with the highest count.
          # rownames and colnames of sdm (all.umi.count) is sorted by read counts
          # increasingly
          sdm.edit.ind <- max(which(sdm[k,] <= umi.edit))
          
          # correct current umi k at position i
          if (which(rownames(sdm) == k) < sdm.edit.ind) {
            rdt[umi == k & position == i, inferred_umi := colnames(sdm)[sdm.edit.ind]]
          }
        }
      }
    }
  }
  return (rdt)
}


na.to.0 = function(dt) {
  # by number (slightly faster than by name)
  for (j in seq_len(ncol(dt)))
    set(dt, which(is.na(dt[[j]])), j, 0)
  return (dt)
}


# correct umi mismatch recursive
umi.mismatch.correction <- function(samdt, current.ref, umi.max.gap, umi.edit) {
  cat(paste0("Processing chr: ", current.ref, "\n"))
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
    cat("Fragment region: \n")
    print(i)
    rdt <- recursive.umi.correction(i, rdt, umi.edit)
  }
  return (rdt)
}


recursive.umi.correction <- function(i, rdt, umi.edit) {
  all.umi.count <- sort(table(c(rdt[position %in% i, inferred_umi])))
  
  if (length(all.umi.count) > 1) {
    sdm <- stringdistmatrix(names(all.umi.count), names(all.umi.count))
    diag(sdm) <- 100
    rownames(sdm) <- names(all.umi.count)
    colnames(sdm) <- names(all.umi.count)
    
    for (j in rownames(sdm)) {
      if (min(sdm[j,]) <= umi.edit) {
        sdm.edit.ind <- max(which(sdm[j,] <= umi.edit))
        # correct current inferred_umi j within position group i
        # if umi j has fewer reads
        if (which(rownames(sdm) == j) < sdm.edit.ind) {
          for (k in rdt[position %in% i & inferred_umi == j, umi]){
            # if edit distance of original umi <= umi.edit
            if (stringdist(k, last(names(which(sdm[j,] <= umi.edit)))) <= umi.edit) {
              rdt[umi == k & inferred_umi == j & position %in% i,
                  inferred_umi := colnames(sdm)[sdm.edit.ind]]
              rdt <- recursive.umi.correction(i, rdt, umi.edit)
              break
            }
          }
        }
      }
    }
  }
  return (rdt)
}


# correct umi mismatch
umi.mismatch.correction.ori <- function(samdt, current.ref, umi.max.gap, umi.edit) {
  # Add inferred_umi info to reference data table
  rdt <- samdt[rname == current.ref,]
  rdt[, c("umi", "inferred_umi") := 
        tstrsplit(qname, ":", fixed=TRUE, 
                  keep=length(tstrsplit(qname, ":", fixed=TRUE)))]
  
  # Correct UMIs with sequencing errors by looking at UMIs in surrounding region
  
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


# INCORRECT number of fragments per transcript function
num.frag.per.transcript <- function(rdt, umi.max.gap) {
  reads.dt <- rdt[,.(rname, inferred_pos, inferred_umi)]
  chr <- rdt$rname[1]
  
  unique.pos <- sort(unique(reads.dt$inferred_pos))
  unique.pos.list <- get.adj.pos.list(unique.pos, umi.max.gap)
  
  res.dt <- data.table(rname = chr, transcript = paste0(chr,"_",1:length(unique.pos.list)))
  
  for (i in 1:length(unique.pos.list)) {
    res.dt[i, num := length(unique(rdt[inferred_pos %in% unique.pos.list[[i]],
                                       inferred_umi]))]
  }
  return (res.dt)
}
