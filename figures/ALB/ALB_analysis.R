setwd('~/Desktop/DriverPower/figures/ALB/')
library(reshape2)
library(ggplot2)
library(ComplexHeatmap)
source('../ggplot_theme.R')
mut = read.table('./Liver-HCC.ALB.mut.tsv', sep='\t', header = F, stringsAsFactors = F)
mut = mut[,c(2,3,4,5,6,8,13)]
colnames(mut) = c('start', 'end', 'type', 'ref', 'alt', 'donor_ID', 'element_ID')
# remove promDomain (all in promCore)
mut = mut[mut$element_ID!='gc19_pc.promDomain::gencode::ALB::ENSG00000163631.12',]
# remove 5utr (all in promCore)
mut = mut[mut$element_ID!='gc19_pc.5utr::gencode::ALB::ENSG00000163631.12',]
table(mut$element_ID)

nrow(unique(mut[,1:6]))
mut = mut[!duplicated(mut[,1:6]),]
mut_simple = dcast(mut, donor_ID ~ element_ID)
colnames(mut_simple) = c('donor_ID', '3utr', 'cds', 'promoter', 'ss')
# CNV
cnv = read.table('./Liver-HCC.ALB.cnv.tsv', sep='\t', header = F, stringsAsFactors = F)
colnames(cnv) = c('donor_ID', 'chrom', 'start', 'end', 'total_cn', 'major_cn', 'minor_cn', 'star')
# for those have CNV break points
part_donors = cnv$donor_ID[duplicated(cnv$donor_ID)]
part_cnv = cnv[cnv$donor_ID %in% part_donors, ]
# for those do not have CNV break points (entire gain of loss)
entire_cnv = cnv[! cnv$donor_ID %in% part_donors, ]

# cnv gain and loss
cnv_simple = entire_cnv[,c(1,5)]
cnv_simple = rbind(cnv_simple,
      c('DO48715', 1.5),
      c('DO48728', NA),
      c('DO45119', 1),
      c('DO50837', 2),
      c('DO50840', 2),
      c('DO45045', 1.5),
      c('DO52171', 2.5),
      c('DO44888', 2),
      c('DO22994', NA),
      c('DO45223', 3.5))
cnv_simple['cnv_gain'] = ifelse(cnv_simple$total_cn>2, 1, 0)
cnv_simple['cnv_loss'] = ifelse(cnv_simple$total_cn<2, 1, 0)

# SV
sv = read.table('./Liver-HCC.ALB.500kb.sv.tsv', header = F, sep='\t', stringsAsFactors = F)
# remove all DUP and DEL events
sv = sv[sv$V12!='DEL',]
sv = sv[sv$V12!='DUP',]
sv_simple = unique(sv[,c(1,12)])
colnames(sv_simple) = c('donor_ID', 'sv_type')
sv_simple['SV'] = 1

# combine
alb = merge(mut_simple, cnv_simple[,c('donor_ID', 'cnv_gain', 'cnv_loss')], all = T)
alb = merge(alb, sv_simple[,c('donor_ID', 'SV')], all=T)
alb[is.na(alb)] = 0
row.names(alb) = alb$donor_ID
alb[alb>0] = 1
alb$donor_ID = row.names(alb)
alb = alb[,c('donor_ID', 'cds', 'ss', '3utr', 'promoter', 'cnv_loss', 'cnv_gain', 'SV')]
alb = alb[order(alb$cds, alb$ss, alb$`3utr`, alb$promoter, alb$cnv_loss, alb$cnv_gain, alb$SV, decreasing = T),]




dat = t(as.matrix(alb[,-1]))
row.names(dat) = c('CDS', 'SS', "3'UTR", 'Promoter', 'CN Loss', 'CN Gain', 'SV')
dat.list = list(dat=dat)
color = c(dat="#fdb462")
png('./ALB.oncoprint.png', width = 960, height = 480)
oncoPrint(dat.list, remove_empty_columns = TRUE,
          alter_fun = list(
            background = function(x, y, w, h) {
              grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.6, "mm"), gp = gpar(fill = "#CCCCCC", col = NA))
            },
            dat = function(x, y, w, h) {
              grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.6, "mm"), gp = gpar(fill = "#386cb0", col = NA))
            }
          ), col=color, heatmap_legend_param = list(), show_row_barplot = F, show_heatmap_legend = FALSE)
dev.off()
write.table(dat, './ALB.oncoprint.txt', sep='\t')

# RNA-seq (187 in total, 69 normals and 118 tumors)
rnameta = read.table('./rnaseq_metadata.lite.tsv', header=T, sep='\t', stringsAsFactors = F)
rnameta = rnameta[rnameta$project_code %in% c('LIHC-US', 'LIRI-JP'),]
table(rnameta$is_tumour)
rnaseq = read.table('~/Downloads/tophat_star_fpkm_uq.v2_aliquot_gl.tsv.gz', header=T, sep='\t', check.names = F)
rnaseq = rnaseq[, c('feature', rnameta$aliquot_id)]
alb.expression = rnaseq[rnaseq$feature=='ENSG00000163631.12',]
alb.expression = t(alb.expression[,-1])
alb.expression = as.data.frame(alb.expression)
alb.expression['aliquot_id'] = row.names(alb.expression)
alb.expression = merge(alb.expression, rnameta[,c(2,3,5)])
colnames(alb.expression)[2] = 'FPKM.UQ'
colnames(alb.expression)[3] = 'donor_ID'
alb.expression$is_tumour = as.factor(alb.expression$is_tumour)
alb.expression = alb.expression[alb.expression$donor_ID %in% alb$donor_ID, ] # some of LIRI-JP are actually Biliary-AdenoCA
write.table(alb.expression, 'ALB.expression.txt', sep='\t')
# expression normal vs. tumor
levels(alb.expression$is_tumour) = c('Normal', 'Tumor')
ylabel = expression(paste("FPKM-UQ expression of ", italic("ALB")))
xlabs <- paste0(levels(alb.expression$is_tumour),"\n(N=",table(alb.expression$is_tumour),")")
ggplot(alb.expression, aes(x=is_tumour, y=FPKM.UQ, fill=is_tumour)) + geom_boxplot() + ylab(ylabel) +
  scale_fill_Publication1() + guides(fill=FALSE) +
  annotate('segment', x=1, xend = 2, y=200000, yend=200000) +
  annotate('text', x=1.5, y=210000, label='p<2.2e-16') +
  theme_Publication() + scale_x_discrete(labels=xlabs) + theme(axis.title.x = element_blank())
ggsave('./ALB.expression.normal.vs.tumor.tiff')
wilcox.test(FPKM.UQ ~ is_tumour, alb.expression, alternative='greater')

# add expression to mutation
alb = merge(alb, alb.expression[alb.expression$is_tumour=='Tumor',c('donor_ID', 'FPKM.UQ')])

cds = cbind(data.frame(rna=alb[alb$cds==1, 9]), name='CDS', stringsAsFactors=F)
ss = cbind(data.frame(rna=alb[alb$ss==1, 9]), name='SS', stringsAsFactors=F)
utr3 = cbind(data.frame(rna=alb[alb$`3utr`==1, 9]), name="3'UTR", stringsAsFactors=F)
prom = cbind(data.frame(rna=alb[alb$promoter==1, 9]), name='Promoter', stringsAsFactors=F)
cnv_loss = cbind(data.frame(rna=alb[alb$cnv_loss==1, 9]), name='CN Loss', stringsAsFactors=F)
cnv_gain = cbind(data.frame(rna=alb[alb$cnv_gain==1, 9]), name='CN Gain', stringsAsFactors=F)
sv = cbind(data.frame(rna=alb[alb$SV==1, 9]), name='SV', stringsAsFactors=F)
wt = cbind(data.frame(rna=alb[rowSums(alb[,2:8])==0, 9]), name='WT', stringsAsFactors=F)
normal = alb.expression[alb.expression$is_tumour=='Normal', c('FPKM.UQ','is_tumour')]
colnames(normal) = c('rna','name')
dat.rna = rbind(cds, ss, utr3, prom, cnv_loss, cnv_gain, sv, wt, normal, stringsAsFactors=F)

# huge boxplot
dat.rna$name = factor(dat.rna$name, levels=c("Normal","WT","CDS","Promoter","SS","3'UTR", 'CN Loss', 'CN Gain', 'SV'))
xlabs <- paste0(levels(dat.rna$name),"\n(N=",table(dat.rna$name),")")
ggplot(dat.rna, aes(x=name, y=rna, fill=name)) + geom_boxplot() + ylab(ylabel) +
  theme_Publication() + scale_x_discrete(labels=xlabs) + theme(axis.title.x = element_blank()) +
  annotate("segment", x=rep(2,3),xend=c(3,7,8), y= c(120000, 140000, 160000), yend=c(120000, 140000, 160000)) +
  annotate("text", x=2.5, y=125000, label="p=0.0175") +
  annotate("text", x=4.5, y=145000, label="p=0.0114") +
  annotate("text", x=5, y=165000, label="p=0.0379") + scale_fill_Publication1() + guides(fill=F)
ggsave('./ALB.expression.wt.vs.mut.tiff')
wilcox.test(wt$rna, cds$rna, alternative='greater') # 0.01746
wilcox.test(wt$rna, prom$rna, alternative='greater') # 0.1044
wilcox.test(wt$rna, utr3$rna, alternative='greater') # 0.2583
wilcox.test(wt$rna, cnv_loss$rna, alternative='greater') # 0.01139
wilcox.test(wt$rna, cnv_gain$rna, alternative='greater') # 0.03791
wilcox.test(wt$rna, sv$rna, alternative='greater') # 0.2389

library(MASS)
test.df = function(x ,y){
  rna = c(x,y)
  grp = c(rep(0, length(x)), rep(1, length(y)))
  anova(glm(rna~grp, family = quasipoisson()), test = 'LRT')
}

# ALB mut maf
maf = read.table('./Liver-HCC.ALB.mut.maf', sep='\t', stringsAsFactors = F, header = F, quote='')
# 675 mutations
# remove mutations in repeat, left 451 mutations
maf = maf[is.na(maf$V11),]
table(maf$V6)
nonsym.donors = unique(maf[maf$V6 %in% c("Frame_Shift_Del", "Frame_Shift_Ins", "In_Frame_Del",
                           "In_Frame_Ins", "Missense_Mutation", "Nonsense_Mutation"), 'V13'])
ss.donors = unique(maf[maf$V6 == 'Splice_Site', 'V13'])
prom.donors = unique(maf[maf$V6 %in% c("5'Flank", "5'UTR"), 'V13'])
