# GPR126
# Bladder-TCC
# No coding mutations
setwd('~/Desktop/DriverPower/figures/')
rnameta = read.table('./ALB/rnaseq_metadata.lite.tsv', header=T, sep='\t', stringsAsFactors = F)
rnameta = rnameta[rnameta$project_code %in% c('BLCA-US'),]
table(rnameta$is_tumour)
rnaseq = read.table('~/Downloads/tophat_star_fpkm_uq.v2_aliquot_gl.tsv.gz', header=T, sep='\t', check.names = F)
rnaseq = rnaseq[, c('feature', rnameta$aliquot_id)]
gpr126.expression = rnaseq[rnaseq$feature=='ENSG00000112414.10',]
gpr126.expression = t(gpr126.expression[,-1])
gpr126.expression = as.data.frame(gpr126.expression)
gpr126.expression['aliquot_id'] = row.names(gpr126.expression)
gpr126.expression = merge(gpr126.expression, rnameta[,c(2,3,5)])
colnames(gpr126.expression)[2] = 'FPKM.UQ'
colnames(gpr126.expression)[3] = 'donor_ID'
gpr126.expression$is_tumour = as.factor(gpr126.expression$is_tumour)
write.table(alb.expression, './GPR126/GPR126.expression.txt', sep='\t')

levels(gpr126.expression$is_tumour) = c('Normal', 'Tumor')
ylabel = expression(paste("FPKM-UQ expression of ", italic("ALB")))
xlabs <- paste0(levels(gpr126.expression$is_tumour),"\n(N=",table(gpr126.expression$is_tumour),")")
ggplot(gpr126.expression, aes(x=is_tumour, y=FPKM.UQ, fill=is_tumour)) + geom_boxplot() + ylab(ylabel) +
  scale_fill_Publication1() + guides(fill=FALSE)
  annotate('segment', x=1, xend = 2, y=200000, yend=200000) +
  annotate('text', x=1.5, y=210000, label='p<2.2e-16') +
  theme_Publication() + scale_x_discrete(labels=xlabs) + theme(axis.title.x = element_blank())
ggsave('./ALB.expression.normal.vs.tumor.tiff')
wilcox.test(FPKM.UQ ~ is_mut, gpr126.expression[gpr126.expression$is_tumour=='Tumor',])

# Donor with mut
donors.in = c('DO44097', "DO472", "DO477", "DO496",
              "DO498", "DO522", "DO555", "DO561",
              "DO669", 'DO695', 'DO804', 'DO822',
              'DO828', 'DO856')
gpr126.expression['is_mut'] = ifelse(gpr126.expression$donor_ID %in% donors.in, 'MUT', 'WT')
gpr126.expression[gpr126.expression$is_tumour=='Normal', 'is_mut'] = 'Normal'
ggplot(gpr126.expression, aes(x=is_mut, y=FPKM.UQ, fill=is_tumour)) + geom_boxplot() + ylab(ylabel) +
  scale_fill_Publication1() + guides(fill=FALSE)
annotate('segment', x=1, xend = 2, y=200000, yend=200000) +
  annotate('text', x=1.5, y=210000, label='p<2.2e-16') +
  theme_Publication() + scale_x_discrete(labels=xlabs) + theme(axis.title.x = element_blank())