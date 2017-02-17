# GPR126
# Bladder-TCC
# No coding mutations
setwd('~/Desktop/DriverPower/figures/')
source('./ggplot_theme.R')
rnameta = read.table('./ALB/rnaseq_metadata.lite.tsv', header=T, sep='\t', stringsAsFactors = F)
rnameta = rnameta[rnameta$project_code %in% c('BLCA-US'),]
table(rnameta$is_tumour)
library(data.table)
library(ggplot2)
rnaseq = fread('~/Downloads/tophat_star_fpkm_uq.v2_aliquot_gl.tsv', select = c('feature', rnameta$aliquot_id))
rnaseq = as.data.frame(rnaseq)
g = 'ENSG00000009844.11' # VTA1
g = 'ENSG00000112414.10' # GPR126
g = 'ENSG00000236366.1' # RP11-440G9.1
g = 'ENSG00000266843.1'
# VEGF pathway genes:
# VEGFR1: gc19_pc.cds::gencode::FLT1::ENSG00000102755.6
# VEGFR2: gc19_pc.cds::gencode::KDR::ENSG00000128052.8
# VEGFR3: gc19_pc.cds::gencode::FLT4::ENSG00000037280.11
# VEGFA: gc19_pc.cds::gencode::VEGFA::ENSG00000112715.16
# VEGFB: gc19_pc.cds::gencode::VEGFB::ENSG00000173511.5
# VEGFC: gc19_pc.cds::gencode::VEGFC::ENSG00000150630.2
# VEGFD: gc19_pc.cds::gencode::FIGF::ENSG00000165197.4
# PGF: gc19_pc.cds::gencode::PGF::ENSG00000119630.9

plot_gene <- function(g, rnaseq, rnameta, donors.in, name){
  g.exp = as.data.frame(t(rnaseq[rnaseq$feature==g,-1]))
  g.exp['aliquot_id'] = row.names(g.exp)
  g.exp = merge(g.exp, rnameta[,c(2,3,5)])
  colnames(g.exp)[2:3] = c('FPKM.UQ', 'donor_ID')
  g.exp['is_mut'] = ifelse(g.exp$donor_ID %in% donors.in, 'MUT', 'WT')
  g.exp[g.exp$is_tumour=='no', 'is_mut'] = 'Normal'
  g.exp$is_mut = factor(g.exp$is_mut, levels = c('Normal', 'WT', 'MUT'))
  xlabs <- paste0(levels(g.exp$is_mut),"\n(N=",table(g.exp$is_mut),")")
  ylabel = paste0("FPKM-UQ expression of ", name)
  ggplot(g.exp, aes(x=is_mut, y=FPKM.UQ, fill=is_tumour)) + geom_boxplot() + ylab(ylabel) +
    scale_fill_Publication1() + guides(fill=FALSE) +
    theme_Publication() + scale_x_discrete(labels=xlabs) + theme(axis.title.x = element_blank())
}
plot_gene('ENSG00000112414.10', rnaseq, rnameta, donors.in, 'GPR126')
# VEGFR
plot_gene('ENSG00000102755.6', rnaseq, rnameta, donors.in, 'FLT1')
plot_gene('ENSG00000128052.8', rnaseq, rnameta, donors.in, 'KDR')
plot_gene('ENSG00000037280.11', rnaseq, rnameta, donors.in, 'FLT4')
# VEGF
plot_gene('ENSG00000112715.16', rnaseq, rnameta, donors.in, 'VEGFA')
plot_gene('ENSG00000173511.5', rnaseq, rnameta, donors.in, 'VEGFB')
plot_gene('ENSG00000150630.2', rnaseq, rnameta, donors.in, 'VEGFC')
plot_gene('ENSG00000165197.4', rnaseq, rnameta, donors.in, 'VEGFD')
plot_gene('ENSG00000119630.9', rnaseq, rnameta, donors.in, 'PGF')



gpr126.expression = rnaseq[rnaseq$feature==g,]
gpr126.expression = t(gpr126.expression[,-1])
gpr126.expression = as.data.frame(gpr126.expression)
gpr126.expression['aliquot_id'] = row.names(gpr126.expression)
gpr126.expression = merge(gpr126.expression, rnameta[,c(2,3,5)])
colnames(gpr126.expression)[2] = 'FPKM.UQ'
colnames(gpr126.expression)[3] = 'donor_ID'
gpr126.expression$is_tumour = as.factor(gpr126.expression$is_tumour)
write.table(gpr126.expression, './GPR126/Bladder-TCC.GPR126.expression.txt', sep='\t')

levels(gpr126.expression$is_tumour) = c('Normal', 'Tumor')
ylabel = expression(paste("FPKM-UQ expression of ", italic("GPR126")))
xlabs <- paste0(levels(gpr126.expression$is_tumour),"\n(N=",table(gpr126.expression$is_tumour),")")
ggplot(gpr126.expression, aes(x=is_tumour, y=FPKM.UQ, fill=is_tumour)) + geom_boxplot() + ylab(ylabel) +
  scale_fill_Publication1() + guides(fill=FALSE) +
  theme_Publication() + scale_x_discrete(labels=xlabs) + theme(axis.title.x = element_blank())
wilcox.test(FPKM.UQ ~ is_mut, gpr126.expression[gpr126.expression$is_mut != 'Normal',])
wilcox.test(FPKM.UQ ~ is_tumour, gpr126.expression)

# Donor with mut
donors.in = c('DO44097', "DO472", "DO477", "DO496",
              "DO498", "DO522", "DO555", "DO561",
              "DO669", 'DO695', 'DO804', 'DO822',
              'DO828', 'DO856')
gpr126.expression['is_mut'] = ifelse(gpr126.expression$donor_ID %in% donors.in, 'MUT', 'WT')
gpr126.expression[gpr126.expression$is_tumour=='Normal', 'is_mut'] = 'Normal'
gpr126.expression$is_mut = factor(gpr126.expression$is_mut, levels = c('Normal', 'WT', 'MUT'))
xlabs <- paste0(levels(gpr126.expression$is_mut),"\n(N=",table(gpr126.expression$is_mut),")")
ggplot(gpr126.expression, aes(x=is_mut, y=FPKM.UQ, fill=is_tumour)) + geom_boxplot() + ylab(ylabel) +
  scale_fill_Publication1() + guides(fill=FALSE) +
  scale_x_discrete(labels=xlabs) + theme(axis.title.x = element_blank())
  annotate('segment', x=2, xend = 3, y=30, yend=30) +
  annotate('text', x=2.5, y=35, label='p=0.04559') +
  theme_Publication()

# Breast Cancer
donor_id = read.table('./GPR126/Breast-AdenoCa.donor_id.txt', stringsAsFactors = F) # 195
rnameta = read.table('./ALB/rnaseq_metadata.lite.tsv', header=T, sep='\t', stringsAsFactors = F)
rnameta = rnameta[rnameta$icgc_donor_id %in% donor_id$V1,]
table(rnameta$is_tumour)
table(rnameta$wgs_exclusion_white_gray)
table(rnameta$wgs_white_black_gray)
rnaseq = fread('~/Downloads/tophat_star_fpkm_uq.v2_aliquot_gl.tsv', select = c('feature', rnameta$aliquot_id))
rnaseq = as.data.frame(rnaseq)
gpr126.expression = rnaseq[rnaseq$feature==g,]
gpr126.expression = t(gpr126.expression[,-1])
gpr126.expression = as.data.frame(gpr126.expression)
gpr126.expression['aliquot_id'] = row.names(gpr126.expression)
gpr126.expression = merge(gpr126.expression, rnameta[,c(2,3,5)])
colnames(gpr126.expression)[2] = 'FPKM.UQ'
colnames(gpr126.expression)[3] = 'donor_ID'
gpr126.expression$is_tumour = as.factor(gpr126.expression$is_tumour)
levels(gpr126.expression$is_tumour) = c('Normal', 'Tumor')
ylabel = expression(paste("FPKM-UQ expression of ", italic("GPR126")))
# Donor with mut
mut = read.table('GPR126/Breast-AdenoCa.gpr126.mut.tsv', stringsAsFactors = F, sep='\t')
donors.in = unique(mut$V8)
gpr126.expression['is_mut'] = ifelse(gpr126.expression$donor_ID %in% donors.in, 'MUT', 'WT')
gpr126.expression[gpr126.expression$is_tumour=='Normal', 'is_mut'] = 'Normal'
gpr126.expression$is_mut = factor(gpr126.expression$is_mut, levels = c('Normal', 'WT', 'MUT'))
xlabs <- paste0(levels(gpr126.expression$is_mut),"\n(N=",table(gpr126.expression$is_mut),")")
ggplot(gpr126.expression, aes(x=is_mut, y=FPKM.UQ, fill=is_tumour)) + geom_boxplot() + ylab(ylabel) +
  scale_fill_Publication1() + guides(fill=FALSE) +
  theme_Publication() + scale_x_discrete(labels=xlabs) + theme(axis.title.x = element_blank())

wilcox.test(FPKM.UQ ~ is_mut, gpr126.expression[gpr126.expression$is_mut!='Normal',])
wilcox.test(FPKM.UQ ~ is_mut, gpr126.expression[gpr126.expression$is_mut!='WT',])
wilcox.test(FPKM.UQ ~ is_mut, gpr126.expression[gpr126.expression$is_mut!='MUT',])


write.table(gpr126.expression, './GPR126/Breast-AdenoCa.GPR126.expression.txt', sep='\t')

# GPR126 exon level expression
rnameta = read.table('./ALB/rnaseq_metadata.lite.tsv', header=T, sep='\t', stringsAsFactors = F)
rnameta = rnameta[rnameta$project_code %in% c('BLCA-US'),]
table(rnameta$is_tumour)
rnameta['mut'] = 'WT'
rnameta[rnameta$icgc_donor_id %in% donors.in, 'mut'] = 'MUT'
rnameta[rnameta$is_tumour=='no', 'mut'] = 'Normal'
rnaseq = read.table('./GPR126/GPR126.exon.expression.tsv', header=T, sep='\t', check.names = F)
rnaseq = rnaseq[, c('feature', rnameta$aliquot_id)]
rnaseq['median.normal'] = apply(rnaseq[, rnameta$aliquot_id[rnameta$is_tumour=='no']], 1, median)
rnaseq['median.mut'] = apply(rnaseq[, rnameta$aliquot_id[rnameta$icgc_donor_id %in% donors.in]], 1, median)
rnaseq['median.wt'] = apply(rnaseq[, rnameta$aliquot_id[rnameta$mut == 'WT']], 1, median)
plot(log2(median.normal+1)~feature, data = rnaseq, las=2)
lines(log2(median.mut+1)~feature, data = rnaseq, col='red')
lines(log2(median.wt+1)~feature, data = rnaseq, col='blue')
