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
process_one <- function(meta, pancan, tumors, name){
  combine = read.table(paste0(tumors[1], '.', name, '.observed.txt'), header=T, stringsAsFactors = F)
  combine = combine[, c('element_ID', 'p.value', 'q.value')]
  colnames(combine) = c('element_ID', paste0('p.', tumors[1]), paste0('q.', tumors[1]))
  combine[is.na(combine)] = 1
  for (tumor in tumors[2:length(tumors)]) {
    tb = read.table(paste0(tumor, '.', name, '.observed.txt'), header=T, stringsAsFactors = F)[, c('element_ID', 'p.value', 'q.value')]
    colnames(tb) = c('element_ID', paste0('p.', tumor), paste0('q.', tumor))
    tb[is.na(tb)] = 1
    combine = merge(combine, tb, by='element_ID')
  }
  combine[, 'N'] = rowSums(combine[,seq(3, 2*length(tumors)+1, 2)]<=0.1) # num of tumors that have this driver
  combine[, 'p.min'] = apply(combine[,seq(2, 2*length(tumors), 2)], 1, min, na.rm=T)
  combine[, 'q.min'] = apply(combine[,seq(3, 2*length(tumors)+1, 2)], 1, min, na.rm=T)
  # read length, nSample
  meta = read.table(meta, header=T, stringsAsFactors = F)
  meta = meta[,c('binID', 'Length', 'nMut', 'nSample')]
  colnames(meta)[1] = 'element_ID'
  combine = merge(combine, meta, by='element_ID', all.x=T)
  # read pancan
  pancan = read.table(pancan, header=T, stringsAsFactors = F)
  pancan = pancan[,c('element_ID', 'p.value', 'q.value')]
  pancan[is.na(pancan)] = 1
  colnames(pancan) = c('element_ID', 'p.PanCan', 'q.PanCan')
  dat = combine[, c('element_ID', 'p.min', 'q.min', 'nSample', 'N')]
  dat = merge(pancan, dat, by='element_ID')
  dat[dat$q.PanCan<1e-20,'q.PanCan'] = 1e-20
  dat[dat$q.min<1e-20,'q.min'] = 1e-20
  dat$x = -log10(dat$q.min)
  dat$y = -log10(dat$q.PanCan)
  return(dat)
}

# cds
setwd('~/Desktop/DriverPower/results/CDS.DriverPower.observed.20170130/')
dat.cds = process_one('../../figures/data/PanCan_No_Melanoma_Lymph.cds.cadd.tsv',
                      './PanCan_No_Melanoma_Lymph.CDS.observed.txt',
                      singleType, 'CDS')
# coding plots
dat.cds['type'] = tstrsplit(dat.cds$element_ID, '::')[[1]]
dat.cds['text'] = tstrsplit(dat.cds$element_ID, '::')[[3]]
sig.only.cds = dat.cds[dat.cds$q.min<=0.1, ]
dat.cds$nSample[is.na(dat.cds$nSample)] = 0
yline.cds = data.frame(y=rep(1,1000), x=seq(0,999))
p.cds = ggplot(data=dat.cds, aes(x=nSample, y=x)) + geom_point(color='grey', size=1) +
  geom_line(data=yline.cds, aes(x, y), color='orange', size=1.5) +
  geom_point(dat=sig.only.cds, aes(x=nSample, y=x, size=N), color="#386cb0") +
  geom_text_repel(dat=sig.only.cds, aes(x=nSample, y=x, label=text), segment.color='blue', size=3, color='black', fontface="italic") +
  coord_trans(x='sqrt', y='sqrt') +
  scale_x_continuous('Number of mutated samples (N=2279)', breaks = seq(0, 901, 100)) + 
  scale_y_continuous('-log10 q value (most significant tumour type)', breaks=seq(1,20,2)) +
  # annotate('text', x=800, y=1.3, label='Sig. hits = 89', size=5) +
  theme_Publication() + scale_color_continuous_tableau() + guides(size=FALSE)
  
ggsave('../../figures/all.cds.qvsn.20170130.tiff', p.cds, height = 10, width = 10, dpi = 300)


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

# promCore
setwd('~/Desktop/DriverPower/results/promCore.DriverPower.observed.20170130//')
dat.prom = process_one('../../figures/data/PanCan_No_Melanoma_Lymph.promCore.cadd.tsv',
                       './PanCan_No_Melanoma_Lymph.promCore.observed.txt',
                       singleType, 'promCore')

# lncrna.ncrna
setwd('~/Desktop/DriverPower/results/lncrna.ncrna.DriverPower.observed.20170130/')
dat.ncrna = process_one('../../figures/data/PanCan_No_Melanoma_Lymph.lncrna.ncrna.cadd.tsv',
                        './PanCan_No_Melanoma_Lymph.lncrna.ncrna.observed.txt',
                        singleType, 'lncrna.ncrna')

# enhancers
setwd('~/Desktop/DriverPower/results/enhancers.DriverPower.observed.20170130/')
dat.enhancers = process_one('../../figures/data/PanCan_No_Melanoma_Lymph.enhancers.cadd.tsv',
                            './PanCan_No_Melanoma_Lymph.enhancers.observed.txt',
                            singleType, 'enhancers')

# 3utr
setwd('~/Desktop/DriverPower/results/3utr.DriverPower.observed.20170130/')
dat.3utr = process_one('../../figures/data/PanCan_No_Melanoma_Lymph.3utr.cadd.tsv',
                       './PanCan_No_Melanoma_Lymph.3utr.observed.txt',
                       singleType, '3utr')

# 5utr
setwd('~/Desktop/DriverPower/results/5utr.DriverPower.observed.20170130/')
dat.5utr = process_one('../../figures/data/PanCan_No_Melanoma_Lymph.5utr.cadd.tsv',
                       './PanCan_No_Melanoma_Lymph.5utr.observed.txt',
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
sig.only[sig.only$element_ID=='enhancers::chr6:142705600-142706400::NA::NA', 'text'] = 'GPR126.enhancer'
sig.only[sig.only$element_ID=='enhancers::chr10:54205800-54213400::NA::NA', 'text'] = 'DKK1/LINC01468.enhancer'
sig.only[sig.only$element_ID=='enhancers::chr7:86865600-86866400::NA::NA', 'text'] = 'TP53TG1.enhancer'
linked.genes = 'ENSG00000107984;ENSG00000231131;ENSG00000107984;ENSG00000107984;ENSG00000107984;ENSG00000107984;ENSG00000107984;ENSG00000107984;ENSG00000107984;ENSG00000107984;ENSG00000107984;ENSG00000231131;ENSG00000107984;ENSG00000107984;ENSG00000107984;ENSG00000107984;ENSG00000107984;ENSG00000107984;ENSG00000231131;ENSG00000231131;ENSG00000107984;ENSG00000107984;ENSG00000107984;ENSG00000107984;ENSG00000107984;ENSG00000231131;ENSG00000107984;ENSG00000107984;ENSG00000107984;ENSG00000107984;ENSG00000107984;ENSG00000107984;ENSG00000107984;ENSG00000107984;ENSG00000107984;ENSG00000231131;ENSG00000231131;ENSG00000231131;ENSG00000231131;ENSG00000231131;ENSG00000107984;ENSG00000231131'
linked.genes = unique(strsplit(linked.genes, ';', fixed = T)[[1]])

dat$nSample[is.na(dat$nSample)] = 0
yline = data.frame(y=rep(1,400), x=seq(0,399))
p.nc = ggplot(data=dat, aes(x=nSample, y=x)) + geom_point(color='grey') +
  geom_line(data=yline, aes(x, y), color='yellow', size=1.5) +
  geom_point(dat=sig.only, aes(x=nSample, y=x, color=type, size=N)) +
  geom_text_repel(dat=sig.only, aes(x=nSample, y=x, label=text), segment.color='black', size=4, color='black',fontface="italic") +
  coord_trans(x='sqrt', y='sqrt') +
  scale_x_continuous('Number of mutated samples (N=2279)') + 
  scale_y_continuous('-log10 q value (most significant tumour type)', breaks=seq(1,20,2)) +
  theme_Publication() + scale_color_Publication1() + guides(color=FALSE, size=FALSE)
ggsave('../../figures/all.non-coding.qvsn.20170130.png', p.nc, height = 10, width = 10, dpi = 300)

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


