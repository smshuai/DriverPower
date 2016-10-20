#!/usr/bin/env Rscript

suppressPackageStartupMessages(require(Biostrings))

args <- commandArgs(trailingOnly = TRUE)

if(length(args) != 3 ) {
    stop("The number of args is incorrect. Usage: fa2cg.R BED FASTA OUTPATH")
    q("no", 1, FALSE)
}

bed     = read.table(args[1])
fa      = readDNAStringSet(args[2], 'fasta')
outpath = args[3]

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

bed = bed[, 1:4]
colnames(bed) = c('chrom', 'start', 'end', 'binID')
output = cbind(bed[, 1:4], coverage)
write.table(output, outpath, sep="\t", row.names=F, quote=F)
