#!/usr/bin/env Rscript

library(stringdist)
library(optparse)

option_list = list(
                                make_option(c("--window"), action="store", type="integer", default=20, help="How many bases to look up/down stream for matching UMIs [default %default]"),
                                make_option(c("--umi_length"), action="store", type="integer", default=8, help="Length of UMI at the end of each read ID [default %default]"),
                                make_option(c("--umi_edit_dist"), action="store", type="integer", default=2, help="Number of mismatches to allow in UMI [default %default]")
                 )
opt = parse_args(OptionParser(usage = "%prog [options] in.sam out.sam", option_list = option_list, description='Parses SAM file (without the header) and marks duplicates with the same Unique Molecular Index (UMI).'),
				print_help_and_exit = TRUE,
				positional_arguments=2)
				
umi.edit = as.numeric(opt$options$umi_edit_dist)
umi.window = as.numeric(opt$options$window)
umi.length = as.numeric(opt$options$umi_length)
in.sam = opt$args[1]
out.sam = opt$args[2]

cat("UMI edit distance:", umi.edit, "\n")
cat("UMI length at the end of the read id:", umi.length, "\n")
cat("Window +/-", umi.window, "bases\n")

cat("\nReading in file", in.sam, "...\n")

duplication.offset = 1024  ## SAM specification for duplication in FLAG column

sam = read.table(in.sam, sep="\t", comment="", quote="", fill=TRUE, header=FALSE, stringsAsFactors=FALSE)
colnames(sam)[1:11] = c("qname", "flag", "rname", "position", "mapq", "cigar", "mref", "mposition", "tlen", "seq", "qual")
chr = setdiff(unique(sam[,3]), "*")

umi.matrix = c()

for(current.ref in chr) {

	cat("Now processing", current.ref, "...\n")

	## Covert to data frame and add UMI/Duplication info
	rdf <- subset(sam, sam$rname == current.ref)

	rdf.umi = substring(rdf$qname, nchar(rdf$qname)-umi.length+1)
	rdf = data.frame(rdf, umi=rdf.umi, inferred_umi=rdf.umi, Duplicate=FALSE, stringsAsFactors=FALSE)

	## Correct UMIs with sequencing errors by looking at UMIs in surrounding region
	unique.pos = sort(unique(rdf$position))
	for(i in unique.pos) {

	  ## Get range data frame for position and surrounding window
	  rdf.sub = subset(rdf, rdf$position == i)
	  rdf.flank.5p = subset(rdf, rdf$position >= (i-umi.window) & rdf$position <= (i-1))
	  rdf.flank.3p = subset(rdf, rdf$position >= (i+1) & rdf$position <= (i+umi.window))
	  rdf.flank = rbind(rdf.flank.5p, rdf.flank.3p)

	  ## Get all unique UMIs in the region
	  all.umi.count = sort(table(c(rdf.sub$umi, rdf.flank$umi)))
	  position.umi.count = sort(table(rdf.sub$umi))

	  ## Align UMIs to all other UMIs. Assign UMIs with lower counts to 
	  ## matching UMI with higher counts
	  if(length(all.umi.count) > 1) {
		sdm = stringdistmatrix(names(all.umi.count), names(all.umi.count))
		diag(sdm) = 100
		rownames(sdm) = names(all.umi.count)
		colnames(sdm) = names(all.umi.count)

		for(k in names(position.umi.count)) {
	
		  ## Get the min edit distance for that UMI
		  sdm.min.edit = min(sdm[k,])
			
		  if(sdm.min.edit <= umi.edit) {
  
			## Get the indices of the edit distances less than umi.edit with the highest count
			sdm.edit.ind = max(which(sdm[k,] <= umi.edit))
			if (which(rownames(sdm) == k) < sdm.edit.ind) {
			  ind = rdf$umi == k & rdf$position == i
			  rdf[ind,"inferred_umi"] = colnames(sdm)[sdm.edit.ind]
			}
		  }
		}
	  }
	}


	## Now that UMIs have been fixed, go through each position and pick one read to
	## represent each UMI. Positions with more reads with that UMI and that are more 5'
	## are given higher priority. One read randomly selected to be the non-duplicate
	for(i in unique.pos) {
	  rdf.sub = subset(rdf, rdf$position == i)
	  rdf.flank.5p = subset(rdf, rdf$position >= (i-umi.window) & rdf$position <= (i-1))
	  rdf.flank.3p = subset(rdf, rdf$position >= (i+1) & rdf$position <= (i+umi.window))
	  rdf.flank = rbind(rdf.flank.5p, rdf.flank.3p)

	  position.umi.count = sort(table(rdf.sub$inferred_umi))

	  for(k in names(position.umi.count)) {
		rdf.umi.sub = subset(rdf.sub, inferred_umi == k)
		rdf.umi.flank.5p = subset(rdf.flank.5p, inferred_umi == k)
		rdf.umi.flank.3p = subset(rdf.flank.3p, inferred_umi == k)

		rdf.umi.flank.5p.table = table(rdf.umi.flank.5p$position)
		rdf.umi.flank.3p.table = table(rdf.umi.flank.3p$position)
	
		rdf.umi.flank.5p.max = ifelse(length(rdf.umi.flank.5p.table) == 0, 0, max(rdf.umi.flank.5p.table))
		rdf.umi.flank.3p.max = ifelse(length(rdf.umi.flank.3p.table) == 0, 0, max(rdf.umi.flank.3p.table))    

		if(nrow(rdf.umi.sub) > rdf.umi.flank.5p.max & nrow(rdf.umi.sub) >= rdf.umi.flank.3p.max) {
		  ind = rdf$qname %in% rdf.umi.sub$qname[-1]
		  rdf[ind,"flag"] = rdf[ind,"flag"] + duplication.offset
		  rdf[ind,"Duplicate"] = TRUE
		} else {
		  ind = rdf$qname %in% rdf.umi.sub$qname
		  rdf[ind,"Duplicate"] = TRUE
		}
	  }
	}
	write.table(rdf[,-match(c("umi", "inferred_umi", "Duplicate"), colnames(rdf))], out.sam, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE)
	
	umi.matrix = rbind(umi.matrix, rdf[,c("rname", "position", "umi", "inferred_umi", "Duplicate")])
}


rdf = subset(sam, rname == '*')
write.table(rdf, out.sam, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE)

stats.label = c("Number of aligned reads:", "Number of duplicates:", "Percent duplication:",
          "Number of unique UMIs:", "Number of unique UMIs after correction:", "Number of reads with mismatches in UMI:",
          "Percentage of reads with mismatches in UMI:")
stats = c(nrow(umi.matrix), sum(umi.matrix$Duplicate), sum(umi.matrix$Duplicate)/nrow(umi.matrix),
		length(unique(umi.matrix[,3])), length(unique(umi.matrix[,4])), sum(umi.matrix[,3] != umi.matrix[,4]),
		sum(umi.matrix[,3] != umi.matrix[,4]) / nrow(umi.matrix))

write.table(cbind(stats.label, stats), "UMI_stats.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

## Calculate relative abundance of UMIs
umi.table = table(umi.matrix[,3])
umi.inferred.table = table(umi.matrix[,4])

## Base percentage of UMIs at each position
base.matrix = c()
for(i in 1:umi.length) {
  base.matrix = cbind(base.matrix, table(substring(unique(umi.matrix[,3]), i, i)))
}
base.matrix.frac = sweep(base.matrix, 2, length(unique(umi.matrix[,3])), "/")

base.inferred.matrix = c()
for(i in 1:umi.length) {
  base.inferred.matrix = cbind(base.inferred.matrix, table(substring(unique(umi.matrix[,4]), i, i)))
}
base.inferred.matrix.frac = sweep(base.inferred.matrix, 2, length(unique(umi.matrix[,4])), "/")


pdf("UMI_QC_plots.pdf", useDingbats=FALSE)
hist(umi.table, freq=FALSE, xlab="Number of reads per UMI", ylab="UMI Fraction", lwd=2, breaks=20, main="Original UMIs")
hist(umi.inferred.table, freq=FALSE, xlab="Number of reads per UMI", ylab="UMI Fraction", lwd=2, breaks=20, main="Inferred UMIs")
barplot(base.matrix.frac, col=1:4, xlab="UMI position", names.arg=(1:umi.length), ylab="Percent")
barplot(base.inferred.matrix.frac, col=1:4, xlab="UMI position", names.arg=(1:umi.length), ylab="Percent")
plot(1:10, type="n")
legend("center", rownames(base.matrix.frac), col=1:4, pch=15)
dev.off()

o = order(umi.table, decreasing=TRUE)
write.table(cbind(names(umi.table), umi.table)[o,], "UMI_original_counts.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
o = order(umi.inferred.table, decreasing=TRUE)
write.table(cbind(names(umi.inferred.table), umi.inferred.table)[o,], "UMI_inferred_counts.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)


umi.loc = apply(umi.matrix[,c(1,2,4)], 1, paste, collapse="_")
umi.loc = gsub(" ", "", umi.loc)
umi.loc.table = table(umi.loc)
o = order(umi.loc.table, decreasing=TRUE)
write.table(cbind(names(umi.loc.table), umi.loc.table)[o,], "UMI_position_counts.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)












