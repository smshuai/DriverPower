library(data.table)
consel = fread('~/Downloads/phastConsElements100way_UCSC_hg19.tsv',sep='\t')
consel[,length:=chromEnd-chromStart]
consel = consel[order(chrom, chromStart), ]
summary(consel$score)
summary(consel$length)
cutoff = median(consel$score) + 2*IQR(consel$score)
chroms = paste0('chr', c(1:22, 'X', 'Y'))
keep = consel[length>=100 & score>cutoff & chrom %in% chroms,]
keep[, 'element_id':=paste0('Cons100VertEl::', chrom, ':', chromStart, '-', chromEnd, '::', score)]
write.table(keep[,.(chrom, chromStart, chromEnd, element_id)], '~/Downloads/Cons100VertEl.bed',
            quote=F, sep='\t', col.names = F, row.names = F)
