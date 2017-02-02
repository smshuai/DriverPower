library(ggplot2)
library(ggrepel)
library(data.table)
library(scales)
source('~/Desktop/DriverPower/figures/ggplot_theme.R')

# No lymph and melanoma
singleType = c(	'Liver-HCC', 'Panc-AdenoCA', 'Prost-AdenoCA', 'Breast-AdenoCa',
                'Kidney-RCC', 'CNS-Medullo', 'Ovary-AdenoCA', 'Eso-AdenoCa',
                'CNS-PiloAstro', 'Panc-Endocrine', 'Stomach-AdenoCA', "Head-SCC",
                "ColoRect-AdenoCA", "Thy-AdenoCA", "Lung-SCC", "Uterus-AdenoCA",
                "Kidney-ChRCC", "Bone-Osteosarc", "CNS-GBM", 'Lung-AdenoCA', "Biliary-AdenoCA",
                'Bone-Leiomyo', 'Bladder-TCC', 'Myeloid-MPN', 'CNS-Oligo', 'Cervix-SCC')
setwd('~/Desktop/DriverPower/results/CDS.DriverPower/')
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
cds.meta = read.table('../../figures/data/PanCan_No_Melanoma_Lymph.cds.cadd.tsv', header=T, stringsAsFactors = F)
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

# Noncoding
# promCore
setwd('~/Desktop/DriverPower/results/promCore.DriverPower/')
combine = read.table('./Liver-HCC.promCore.DriverPower.observed.txt', header=T, stringsAsFactors = F)
combine = combine[, 4:6]
colnames(combine) = c('element_ID', 'p.Liver-HCC', 'q.Liver-HCC')
for (tumor in singleType[2:26]) {
  tb = read.table(paste0(tumor, '.promCore.DriverPower.observed.txt'), header=T, stringsAsFactors = F)[, 4:6]
  colnames(tb) = c('element_ID', paste0('p.', tumor), paste0('q.', tumor))
  combine = merge(combine, tb, by='element_ID')
}
combine[, 'p.min'] = apply(combine[,seq(2, 52, 2)], 1, min)
combine[, 'q.min'] = apply(combine[,seq(3, 53, 2)], 1, min)
# read length, nSample
meta = read.table('../../figures/data/PanCan_No_Melanoma_Lymph.promCore.cadd.tsv', header=T, stringsAsFactors = F)
meta = meta[,c('binID', 'Length', 'nMut', 'nSample')]
colnames(meta)[1] = 'element_ID'
combine = merge(combine, meta, by='element_ID', all.x=T)

# read pancan
pancan = read.table('./PanCan_No_Melanoma_Lymph.promCore.DriverPower.observed.txt', header=T, stringsAsFactors = F)
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
dat$element_ID = paste0(dat$element_ID, '.p')
sig.only = dat[dat$q.min<=0.1 | dat$q.PanCan <= 0.1, ]

p1 = ggplot(data=dat, aes(x=x, y=y)) + geom_point(color=dat$pcolor) +
  geom_line(data=xline, aes(x, y), color='yellow', size=1.5) + 
  geom_line(data=yline, aes(x, y), color='yellow', size=1.5) +
  geom_text_repel(dat=sig.only, aes(x=x, y=y, label=element_ID), color=sig.only$pcolor, fontface="italic") +
  coord_trans(x='sqrt', y='sqrt') +
  scale_x_continuous('-log10 q value (most significant tumour type)', breaks=seq(1,20,2)) + 
  scale_y_continuous('-log10 q value (pan-cancer cohort)', breaks=seq(1,20,2)) +
  theme_Publication()
ggsave('../../figures/all.promCore.scatter.png', p1, width = 8, height = 8, dpi=500)

dat$nSample[is.na(dat$nSample)] = 0
yline = data.frame(y=rep(1,100), x=seq(0,99))
sig.only = dat[dat$q.min<=0.1, ]
p2 = ggplot(data=dat, aes(x=nSample, y=x)) + geom_point(color=dat$pcolor) +
  geom_line(data=yline, aes(x, y), color='yellow', size=1.5) +
  geom_text_repel(dat=sig.only, aes(x=nSample, y=x, label=element_ID), color=sig.only$pcolor, fontface="italic") +
  coord_trans(x='identity', y='sqrt') +
  scale_x_continuous('Number of mutated samples (N=2279)') + 
  scale_y_continuous('-log10 q value (most significant tumour type)', breaks=seq(1,20,2)) +
  theme_Publication()
ggsave('../../figures/all.promCore.scatter.png', p1, width = 8, height = 8, dpi=500)

# ncrna
setwd('~/Desktop/DriverPower/results/lncrna.ncrna.DriverPower/')
combine.ncrna = read.table('./Liver-HCC.lncrna.ncrna.DriverPower.observed.txt', header=T, stringsAsFactors = F)
combine.ncrna = combine.ncrna[, 4:6]
colnames(combine.ncrna) = c('element_ID', 'p.Liver-HCC', 'q.Liver-HCC')
for (tumor in singleType[2:26]) {
  tb = read.table(paste0(tumor, '.lncrna.ncrna.DriverPower.observed.txt'), header=T, stringsAsFactors = F)[, 4:6]
  colnames(tb) = c('element_ID', paste0('p.', tumor), paste0('q.', tumor))
  combine.ncrna = merge(combine.ncrna, tb, by='element_ID')
}
combine.ncrna[, 'p.min'] = apply(combine.ncrna[,seq(2, 52, 2)], 1, min)
combine.ncrna[, 'q.min'] = apply(combine.ncrna[,seq(3, 53, 2)], 1, min)

# read length, nSample
meta = read.table('../../figures/data/PanCan_No_Melanoma_Lymph.lncrna.ncrna.cadd.tsv', header=T, stringsAsFactors = F)
meta = meta[,c('binID', 'Length', 'nMut', 'nSample')]
colnames(meta)[1] = 'element_ID'
combine.ncrna = merge(combine.ncrna, meta, by='element_ID', all.x=T)

# read pancan
pancan = read.table('./PanCan_No_Melanoma_Lymph.lncrna.ncrna.DriverPower.observed.txt', header=T, stringsAsFactors = F)
pancan = pancan[,4:6]
colnames(pancan) = c('element_ID', 'p.PanCan', 'q.PanCan')
dat.ncrna = combine.ncrna[, c('element_ID', 'p.min', 'q.min', 'nSample')]
dat.ncrna = merge(pancan, dat.ncrna, by='element_ID')
dat.ncrna[dat.ncrna$q.PanCan<1e-20,'q.PanCan'] = 1e-20
dat.ncrna[dat.ncrna$q.min<1e-20,'q.min'] = 1e-20
dat.ncrna$x = -log10(dat.ncrna$q.min)
dat.ncrna$y = -log10(dat.ncrna$q.PanCan)
# color
dat.ncrna['pcolor'] = 'grey'
dat.ncrna[dat.ncrna$q.PanCan<=0.1, 'pcolor'] = 'red'
dat.ncrna[dat.ncrna$q.min<=0.1, 'pcolor'] = 'blue'
dat.ncrna[dat.ncrna$q.min<=0.1 & dat.ncrna$q.PanCan<=0.1, 'pcolor'] = 'black'
xline = data.frame(x=rep(1,22), y=seq(0,21))
yline = data.frame(y=rep(1,22), x=seq(0,21))
dat.ncrna$element_ID = tstrsplit(dat.ncrna$element_ID, split='::', fix=T)[[3]]

# combine ncrna and prom
dat = rbind(dat, dat.ncrna)


# 3utr
setwd('~/Desktop/DriverPower/results/3utr.DriverPower/')
combine.3utr = read.table('./Liver-HCC.3utr.DriverPower.observed.txt', header=T, stringsAsFactors = F)
combine.3utr = combine.3utr[, 4:6]
colnames(combine.3utr) = c('element_ID', 'p.Liver-HCC', 'q.Liver-HCC')
for (tumor in singleType[2:26]) {
  tb = read.table(paste0(tumor, '.3utr.DriverPower.observed.txt'), header=T, stringsAsFactors = F)[, 4:6]
  colnames(tb) = c('element_ID', paste0('p.', tumor), paste0('q.', tumor))
  combine.3utr = merge(combine.3utr, tb, by='element_ID')
}
combine.3utr[, 'p.min'] = apply(combine.3utr[,seq(2, 52, 2)], 1, min)
combine.3utr[, 'q.min'] = apply(combine.3utr[,seq(3, 53, 2)], 1, min)

# read length, nSample
meta = read.table('../../figures/data/PanCan_No_Melanoma_Lymph.3utr.cadd.tsv', header=T, stringsAsFactors = F)
meta = meta[,c('binID', 'Length', 'nMut', 'nSample')]
colnames(meta)[1] = 'element_ID'
combine.3utr = merge(combine.3utr, meta, by='element_ID', all.x=T)

# read pancan
pancan = read.table('./PanCan_No_Melanoma_Lymph.3utr.DriverPower.observed.txt', header=T, stringsAsFactors = F)
pancan = pancan[,4:6]
colnames(pancan) = c('element_ID', 'p.PanCan', 'q.PanCan')
dat.3utr = combine.3utr[, c('element_ID', 'p.min', 'q.min', 'nSample')]
dat.3utr = merge(pancan, dat.3utr, by='element_ID')
dat.3utr[dat.3utr$q.PanCan<1e-20,'q.PanCan'] = 1e-20
dat.3utr[dat.3utr$q.min<1e-20,'q.min'] = 1e-20
dat.3utr$x = -log10(dat.3utr$q.min)
dat.3utr$y = -log10(dat.3utr$q.PanCan)
# color
dat.ncrna['pcolor'] = 'grey'
dat.ncrna[dat.ncrna$q.PanCan<=0.1, 'pcolor'] = 'red'
dat.ncrna[dat.ncrna$q.min<=0.1, 'pcolor'] = 'blue'
dat.ncrna[dat.ncrna$q.min<=0.1 & dat.ncrna$q.PanCan<=0.1, 'pcolor'] = 'black'
xline = data.frame(x=rep(1,22), y=seq(0,21))
yline = data.frame(y=rep(1,22), x=seq(0,21))
dat.ncrna$element_ID = tstrsplit(dat.ncrna$element_ID, split='::', fix=T)[[3]]

# combine ncrna and prom
dat = rbind(dat, dat.ncrna)

process_one <- function(meta, pancan, tumors, name){
  combine = read.table(paste0(tumors[1], '.', name, '.DriverPower.observed.txt'), header=T, stringsAsFactors = F)
  combine = combine[, 4:6]
  colnames(combine) = c('element_ID', paste0('p.', tumors[1]), paste0('q.', tumors[1]))
  for (tumor in tumors[2:length(tumors)]) {
    tb = read.table(paste0(tumor, '.', name, '.DriverPower.observed.txt'), header=T, stringsAsFactors = F)[, 4:6]
    colnames(tb) = c('element_ID', paste0('p.', tumor), paste0('q.', tumor))
    combine = merge(combine, tb, by='element_ID')
  }
  combine[, 'p.min'] = apply(combine[,seq(2, 2*length(tumors), 2)], 1, min)
  combine[, 'q.min'] = apply(combine[,seq(3, 2*length(tumors)+1, 2)], 1, min)
  # read length, nSample
  meta = read.table(meta, header=T, stringsAsFactors = F)
  meta = meta[,c('binID', 'Length', 'nMut', 'nSample')]
  colnames(meta)[1] = 'element_ID'
  combine = merge(combine, meta, by='element_ID', all.x=T)
  # read pancan
  pancan = read.table(pancan, header=T, stringsAsFactors = F)
  pancan = pancan[,4:6]
  colnames(pancan) = c('element_ID', 'p.PanCan', 'q.PanCan')
  dat = combine[, c('element_ID', 'p.min', 'q.min', 'nSample')]
  dat = merge(pancan, dat, by='element_ID')
  dat[dat$q.PanCan<1e-20,'q.PanCan'] = 1e-20
  dat[dat$q.min<1e-20,'q.min'] = 1e-20
  dat$x = -log10(dat$q.min)
  dat$y = -log10(dat$q.PanCan)
  return(dat)
}

# cds
setwd('~/Desktop/DriverPower/results/CDS.DriverPower/')
dat.cds = process_one('../../figures/data/PanCan_No_Melanoma_Lymph.cds.cadd.tsv',
                      './PanCan_No_Melanoma_Lymph.cds.DriverPower.observed.txt',
                      singleType, 'cds')

# promCore
setwd('~/Desktop/DriverPower/results/promCore.DriverPower/')
dat.prom = process_one('../../figures/data/PanCan_No_Melanoma_Lymph.promCore.cadd.tsv',
                       './PanCan_No_Melanoma_Lymph.promCore.DriverPower.observed.txt',
                       singleType, 'promCore')

# lncrna.ncrna
setwd('~/Desktop/DriverPower/results/lncrna.ncrna.DriverPower/')
dat.ncrna = process_one('../../figures/data/PanCan_No_Melanoma_Lymph.lncrna.ncrna.cadd.tsv',
                        './PanCan_No_Melanoma_Lymph.lncrna.ncrna.DriverPower.observed.txt',
                        singleType, 'lncrna.ncrna')

# enhancers
setwd('~/Desktop/DriverPower/results/enhancers.DriverPower/')
dat.enhancers = process_one('../../figures/data/PanCan_No_Melanoma_Lymph.enhancers.cadd.tsv',
                            './PanCan_No_Melanoma_Lymph.enhancers.DriverPower.observed.txt',
                            singleType, 'enhancers')

# 3utr
setwd('~/Desktop/DriverPower/results/3utr.DriverPower/')
dat.3utr = process_one('../../figures/data/PanCan_No_Melanoma_Lymph.3utr.cadd.tsv',
                       './PanCan_No_Melanoma_Lymph.3utr.DriverPower.observed.txt',
                       singleType, '3utr')

# 5utr
setwd('~/Desktop/DriverPower/results/5utr.DriverPower/')
dat.5utr = process_one('../../figures/data/PanCan_No_Melanoma_Lymph.5utr.cadd.tsv',
                       './PanCan_No_Melanoma_Lymph.5utr.DriverPower.observed.txt',
                       singleType, '5utr')

# noncoding plots
dat = rbind(dat.prom, dat.ncrna, dat.3utr, dat.5utr, dat.enhancers)
dat['type'] = tstrsplit(dat$element_ID, '::')[[1]]
dat['text'] = tstrsplit(dat$element_ID, '::')[[3]]
dat['text'] = apply(dat, 1, function(x) 
  if(x['type'] == 'gc19_pc.promCore')
  {
    paste0(x['text'], '.prom')
  }
  else if(x['type'] == "gc19_pc.3utr")
  {
    paste0(x['text'], '.3utr')
  }
  else if(x['type'] == 'gc19_pc.5utr')
  {
    paste0(x['text'], '.5utr')  
  }
  else
  {
    x['text']
  }
)
sig.only = dat[dat$q.min<=0.1, ]
sig.only[sig.only$type=='enhancers', 'text'] = 'ADGRG6.enhancer'
dat$nSample[is.na(dat$nSample)] = 0
yline = data.frame(y=rep(1,400), x=seq(0,399))
p.nc = ggplot(data=dat, aes(x=nSample, y=x)) + geom_point(color='grey') +
  geom_line(data=yline, aes(x, y), color='yellow', size=1.5) +
  geom_point(dat=sig.only, aes(x=nSample, y=x, color=type)) +
  geom_text_repel(dat=sig.only, aes(x=nSample, y=x, label=text, color=type), fontface="italic") +
  coord_trans(x='sqrt', y='sqrt') +
  scale_x_continuous('Number of mutated samples (N=2279)') + 
  scale_y_continuous('-log10 q value (most significant tumour type)', breaks=seq(1,20,2)) +
  theme_Publication() + scale_color_Publication1() + guides(color=FALSE)
ggsave('../../figures/all.non-coding.qvsn.png', p.nc, height = 8, width = 8, dpi = 600)

xline = data.frame(x=rep(1,22), y=seq(0,21))
sig.only = dat[dat$q.min<=0.1 | dat$q.PanCan <= 0.1, ]
sig.only[sig.only$q.PanCan<=0.1, 'pcolor'] = 'pancan'
sig.only[sig.only$q.min<=0.1, 'pcolor'] = 'single'
sig.only[sig.only$q.min<=0.1 & sig.only$q.PanCan<=0.1, 'pcolor'] = 'both'
sig.only[sig.only$element_ID=='enhancers::chr6:142705600-142706400::NA::NA', 'text'] = 'ADGRG6.enhancer'
sig.only[sig.only$element_ID=='enhancers::chr7:86865600-86866400::NA::NA', 'text'] = 'TP53TG1.enhancer'
p.nc2 = ggplot(data=dat, aes(x=y, y=x)) + geom_point(color='grey') +
  geom_line(data=xline, aes(x, y), color='yellow', size=1.5) + 
  geom_line(data=xline, aes(y, x), color='yellow', size=1.5) +
  geom_point(dat=sig.only, aes(x=y, y=x, color=pcolor)) +
  geom_text_repel(dat=sig.only, aes(x=y, y=x, label=text, color=pcolor), fontface="italic") +
  coord_trans(x='sqrt', y='sqrt') +
  scale_x_continuous('-log10 q value (pan-cancer cohort)', breaks=seq(1,20,2)) + 
  scale_y_continuous('-log10 q value (most significant tumour type)', breaks=seq(1,20,2)) +
  theme_Publication()  +  guides(color=FALSE)
ggsave('../../figures/all.non-coding.scatter.png', p.nc2, width = 8, height = 8, dpi=600)


# coding plots
dat.cds['type'] = tstrsplit(dat.cds$element_ID, '::')[[1]]
dat.cds['text'] = tstrsplit(dat.cds$element_ID, '::')[[3]]
sig.only.cds = dat.cds[dat.cds$q.min<=0.1, ]
dat.cds$nSample[is.na(dat.cds$nSample)] = 0
yline.cds = data.frame(y=rep(1,1000), x=seq(0,999))
p.cds = ggplot(data=dat.cds, aes(x=nSample, y=x)) + geom_point(color='grey') +
  geom_line(data=yline.cds, aes(x, y), color='yellow', size=1.5) +
  geom_point(dat=sig.only.cds, aes(x=nSample, y=x, color=type)) +
  geom_text_repel(dat=sig.only.cds, aes(x=nSample, y=x, label=text, color=type), fontface="italic") +
  coord_trans(x='sqrt', y='sqrt') +
  scale_x_continuous('Number of mutated samples (N=2279)', breaks = seq(0, 901, 100)) + 
  scale_y_continuous('-log10 q value (most significant tumour type)', breaks=seq(1,20,2)) +
  theme_Publication() + scale_color_Publication1() + guides(color=FALSE)
ggsave('../../figures/all.cds.qvsn.png', p.cds, height = 8, width = 8, dpi = 600)


xline = data.frame(x=rep(1,22), y=seq(0,21))
sig.only.cds = dat.cds[dat.cds$q.min<=0.1 | dat.cds$q.PanCan <= 0.1, ]
sig.only.cds[sig.only.cds$q.PanCan<=0.1, 'pcolor'] = 'pancan'
sig.only.cds[sig.only.cds$q.min<=0.1, 'pcolor'] = 'single'
sig.only.cds[sig.only.cds$q.min<=0.1 & sig.only.cds$q.PanCan<=0.1, 'pcolor'] = 'both'
p.cds2 = ggplot(data=dat.cds, aes(x=y, y=x)) + geom_point(color='grey') +
  geom_line(data=xline, aes(x, y), color='yellow', size=1.5) + 
  geom_line(data=xline, aes(y, x), color='yellow', size=1.5) +
  geom_point(dat=sig.only.cds, aes(x=y, y=x, color=pcolor)) +
  geom_text_repel(dat=sig.only.cds, aes(x=y, y=x, label=text, color=pcolor), fontface="italic") +
  coord_trans(x='sqrt', y='sqrt') +
  scale_x_continuous('-log10 q value (pan-cancer cohort)', breaks=seq(1,20,2)) + 
  scale_y_continuous('-log10 q value (most significant tumour type)', breaks=seq(1,20,2)) +
  theme_Publication()  + guides(color=FALSE)
ggsave('../../figures/all.cds.scatter.png', p.cds2, width = 8, height = 8, dpi=600)
