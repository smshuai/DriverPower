library(ggplot2)
library(ggrepel)
library(data.table)
library(scales)
setwd('~/DriverPower/results/CDS.DriverPower/')
source('../../figures/ggplot_theme.R')
# No lymph and melanoma
singleType = c(	'Liver-HCC', 'Panc-AdenoCA', 'Prost-AdenoCA', 'Breast-AdenoCa',
                'Kidney-RCC', 'CNS-Medullo', 'Ovary-AdenoCA', 'Eso-AdenoCa',
                'CNS-PiloAstro', 'Panc-Endocrine', 'Stomach-AdenoCA', "Head-SCC",
                "ColoRect-AdenoCA", "Thy-AdenoCA", "Lung-SCC", "Uterus-AdenoCA",
                "Kidney-ChRCC", "Bone-Osteosarc", "CNS-GBM", 'Lung-AdenoCA', "Biliary-AdenoCA",
                'Bone-Leiomyo', 'Bladder-TCC', 'Myeloid-MPN', 'CNS-Oligo', 'Cervix-SCC')
combine = read.table('./Liver-HCC.cds.DriverPower.observed.txt', header=T, stringsAsFactors = F)
combine = combine[, 4:6]
colnames(combine) = c('element_ID', 'p.Liver-HCC', 'q.Liver-HCC')
for (tumor in singleType[2:26]) {
  tb = read.table(paste0(tumor, '.cds.DriverPower.observed.txt'), header=T, stringsAsFactors = F)[, 4:6]
  colnames(tb) = c('element_ID', paste0('p.', tumor), paste0('q.', tumor))
  combine = merge(combine, tb, by='element_ID')
}
combine[, 'p.min'] = apply(combine[,seq(2, 52, 2)], 1, min)
combine[, 'q.min'] = apply(combine[,seq(3, 53, 2)], 1, min)

# read length, nMut and nSample
cds.meta = read.table('../../figures/PanCan_No_Melanoma_Lymph.cds.cadd.tsv', header=T, stringsAsFactors = F)
cds.meta = cds.meta[,c('binID', 'Length', 'nMut', 'nSample')]
colnames(cds.meta)[1] = 'element_ID'
combine = merge(combine, cds.meta, by='element_ID', all.x=T)

dat = combine[, c('element_ID', 'p.min', 'q.min', 'nSample')]
dat$element_ID = tstrsplit(dat$element_ID, split='::', fix=T)[[3]]
dat = dat[order(dat$p.min), ]

ggplot(dat, aes(x=nSample/2279, y=neglog10q)) + geom_point() + geom_abline(intercept = 1, slope = 0)

# read pancan
pancan = read.table('./PanCan_No_Melanoma_Lymph.cds.DriverPower.observed.txt', header=T, stringsAsFactors = F)
pancan = pancan[,4:6]
colnames(pancan) = c('element_ID', 'p.PanCan', 'q.PanCan')
dat = combine[, c('element_ID', 'p.min', 'q.min', 'nSample')]
dat = merge(pancan, dat, by='element_ID')
dat[dat$q.PanCan<1e-20,'q.PanCan'] = 1e-20
dat[dat$q.min<1e-20,'q.min'] = 1e-20
dat$x = -log10(dat$q.min)
dat$y = -log10(dat$q.PanCan)
# color
dat['pcolor'] = 'grey'
dat[dat$q.PanCan<=0.1, 'pcolor'] = 'red'
dat[dat$q.min<=0.1, 'pcolor'] = 'blue'
dat[dat$q.min<=0.1 & dat$q.PanCan<=0.1, 'pcolor'] = 'black'
xline = data.frame(x=rep(1,22), y=seq(0,21))
yline = data.frame(y=rep(1,22), x=seq(0,21))
dat$element_ID = tstrsplit(dat$element_ID, split='::', fix=T)[[3]]
sig.only = dat[dat$q.min<=0.1 | dat$q.PanCan <= 0.1, ]

p = ggplot(data=dat, aes(x=x, y=y)) + geom_point(color=dat$pcolor) +
  geom_line(data=xline, aes(x, y), color='yellow') + 
  geom_line(data=yline, aes(x, y), color='yellow') +
  geom_text_repel(dat=sig.only, aes(x=x, y=y, label=element_ID), color=sig.only$pcolor, fontface="italic") +
  coord_trans(x='sqrt', y='sqrt') +
  scale_x_continuous('-log10 q value (most significant tumour type)', breaks=seq(1,20,2)) + 
  scale_y_continuous('-log10 q value (pan-cancer cohort)', breaks=seq(1,20,2)) +
  theme_Publication()
ggsave('../../figures/all.cds.scatter.png', p, width = 8, height = 8, dpi=500)
