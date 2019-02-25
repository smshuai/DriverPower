#!/usr/bin/env Rscript

suppressPackageStartupMessages(require(Biostrings))
suppressPackageStartupMessages(require(data.table))
args <- commandArgs(trailingOnly = TRUE)

if(length(args) != 5 ) {
    stop("The number of args is incorrect. Usage: fa2cg.R BED BEDraw FASTA CGPATH TOTCGPATH")
    q("no", 1, FALSE)
}

bed     = read.table(args[1])
bedraw  = read.table(args[2])
fa      = readDNAStringSet(args[3], 'fasta')
cgpath    = args[4]
totcgpath = args[5]

bed = bed[, 1:4]
colnames(bed) = c('chrom', 'start', 'end', 'binID')
if (identical(names(fa), paste0(bed$chrom, ':', bed$start, '-', bed$end))) {
    cat('INFO | Sequence names checked\n')
} else {
    stop("Sequence names in FASTA do not match BED. Aborting")
    q("no", 1, FALSE)
}

# 64 sequence contexts
tri_df = trinucleotideFrequency(fa)

# 32 sequence contexts
context = c("ACA", "ACC", "ACG", "ACT", "ATA", "ATC", "ATG", "ATT", "CCA", "CCC", "CCG", "CCT", "CTA", "CTC", "CTG", "CTT", "GCA", "GCC", "GCG", "GCT", "GTA", "GTC", "GTG", "GTT", "TCA", "TCC", "TCG", "TCT", "TTA", "TTC", "TTG", "TTT")
re_context = sapply(context, function(x){as.character(reverseComplement(DNAString(x)))}) 

coverage = matrix(nrow=nrow(tri_df), ncol=32)
colnames(coverage) = context
for (i in 1:ncol(coverage)) {
	tri_nu    = context[i] # 3-mer
	re_tri_nu = re_context[tri_nu] # reverse complement 3-mer
	coverage[,i] = tri_df[, tri_nu] + tri_df[, re_tri_nu] # sum
}

coverage = as.data.table(coverage)
coverage[, binID:=bed$binID]
coverage = coverage[, lapply(.SD, sum), by=binID]
setkey(coverage, binID)
# missing regions
uniq.bins = unique(as.character(bedraw$V4))
nmiss = length(uniq.bins) - nrow(coverage)
if (nmiss > 0) {
    cat('INFO | Fill', nmiss, 'NA bins with 0\n')
    coverage = coverage[uniq.bins]
    for (i in seq_along(coverage)) set(coverage, i=which(is.na(coverage[[i]])), j=i, value=0)
}
setorder(coverage, binID)
totcg = data.frame(binID=coverage$binID, totcg=rowSums(coverage[,2:33, with=F]))


write.table(coverage, cgpath, sep="\t", row.names=F, quote=F)
write.table(totcg, totcgpath, sep='\t', row.names=F, quote=F)
