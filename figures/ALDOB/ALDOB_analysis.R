setwd('~/Desktop/DriverPower/figures/ALDOB/')
# ALDOB mut maf
maf = read.table('./Liver-HCC.ALDOB.mut.maf', sep='\t', stringsAsFactors = F, header = F, quote='')
# 100 mutations
# remove mutations in repeat, left 70 mutations
maf = maf[maf$V33=='',]
table(maf$V6)
nonsym.donors = unique(maf[maf$V6 %in% c("Frame_Shift_Del", "Frame_Shift_Ins", "In_Frame_Del",
                                         "In_Frame_Ins", "Missense_Mutation", "Nonsense_Mutation", 'Splice_Site'), c('V43','V6')])
intron.donors = unique(maf[maf$V6 %in% c("Intron"), c('V43','V6')])
# CNV
cnv = read.table('./Liver-HCC.ALDOB.cnv.tsv', sep='\t', header = F, stringsAsFactors = F)
colnames(cnv) = c('donor_ID', 'chrom', 'start', 'end', 'total_cn', 'major_cn', 'minor_cn', 'star')
cnv = cnv[!duplicated(cnv$donor_ID),] # "DO44890" has two records, both with CN=3
cnv['cn_gain'] = ifelse(cnv$total_cn>2, 1, 0) # 71 donors
cnv['cn_loss'] = ifelse(cnv$total_cn<2, 1, 0) # 42 donors

# merge
aldob = cnv[,c('donor_ID', 'cn_gain', 'cn_loss')]
aldob['nonsyn'] = ifelse(aldob$donor_ID %in% nonsym.donors$V43, 1, 0)
aldob['intron'] = ifelse(aldob$donor_ID %in% intron.donors$V43, 1, 0)
row.names(aldob) = aldob$donor_ID
dat = t(as.matrix(aldob[,-1]))
row.names(dat) = c('CN Gain', 'CN Loss', "Nonsyn. SSM", 'Intron SSM')
dat.list = list(dat=dat)
color = c(dat="#fdb462")
png('./ALDOB.oncoprint.png', width = 960, height = 300)
oncoPrint(dat.list, remove_empty_columns = TRUE,
          alter_fun = list(
            background = function(x, y, w, h) {
              grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#CCCCCC", col = NA))
            },
            dat = function(x, y, w, h) {
              grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#386cb0", col = NA))
            }
          ), col=color, heatmap_legend_param = list(), show_row_barplot = F, show_heatmap_legend = FALSE, row_order = NULL)
dev.off()
write.table(dat, './ALDOB.oncoprint.txt', sep='\t')

# RNA-seq
library(data.table)
rnameta = read.table('../ALB/rnaseq_metadata.lite.tsv', header=T, sep='\t', stringsAsFactors = F)
donor_id = read.table('../ALB/Liver-HCC.tumor_barcode.txt', header=F, sep=' ', stringsAsFactors = F)
rnameta = rnameta[rnameta$icgc_donor_id %in% donor_id$V2,]
pcawg = fread('~/Downloads/tophat_star_fpkm_uq.v2_aliquot_gl.tsv', select = c('feature', rnameta$aliquot_id))
pcawg = as.data.frame(pcawg)
g.expression = pcawg[pcawg$feature=='ENSG00000136872.13',]
g.expression = as.data.frame(t(g.expression[,-1]))
g.expression['aliquot_id'] = row.names(g.expression)
g.expression = merge(g.expression, rnameta[,c(2,3,5)])
colnames(g.expression)[2:3] = c('FPKM.UQ', 'donor_ID')
# ADD mutation type
g.expression['mut'] = 'WT'
g.expression[g.expression$donor_ID %in% nonsym.donors$V43, 'mut'] = 'CDS'
# CNV
cn_gain.donors = aldob$donor_ID[aldob$cn_gain > 0]
g.expression[g.expression$donor_ID %in% cn_gain.donors, 'mut'] = 'CN Gain'
cn_loss.donors = aldob$donor_ID[aldob$cn_loss > 0]
g.expression[g.expression$donor_ID %in% cn_loss.donors, 'mut'] = 'CN Loss'
# normal
g.expression[g.expression$is_tumour=='no', 'mut'] = 'Normal'

# plot
ylabel = expression(paste("FPKM-UQ expression of ", italic("ALDOB")))
g.expression$mut = factor(g.expression$mut, levels = c('Normal', 'WT', 'CDS', 'CN Gain', 'CN Loss'))
xlabs <- paste0(levels(g.expression$mut),"\n(N=",table(g.expression$mut),")")
ggplot(g.expression, aes(x=mut, y=FPKM.UQ, fill=mut)) + geom_boxplot() + ylab(ylabel) +
  theme_Publication() + scale_x_discrete(labels=xlabs) + theme(axis.title.x = element_blank()) +
  scale_fill_Publication1() + guides(fill=F)
wilcox.test(FPKM.UQ ~ mut, data = g.expression[g.expression$mut %in% c('CN Gain', 'WT'),])
wilcox.test(FPKM.UQ ~ mut, data = g.expression[g.expression$mut %in% c('CN Loss', 'WT'),])
wilcox.test(FPKM.UQ ~ mut, data = g.expression[g.expression$mut %in% c('CDS', 'WT'),])
