setwd('~/DriverPower/')
library(ggplot2)
library(data.table)
source('figures/ggplot_theme.R')
meta = read.table('./results/feature_meta.v2.tsv', header=T, sep='\t')
meta = meta[, c('feature_name', 'description')]
meta['short'] = tstrsplit(meta$description, ' | ', fix=T)[[1]]
# ICGC_ESAD
esad = read.table('./figures/data/ICGC_ESAD.select.tsv', header=T, sep='\t')
esad = merge(esad, meta, by.x='fname', by.y='feature_name')
esad = esad[order(esad$rndlasso, decreasing = T),]
esad.top10 = esad[1:10,]
ggplot(esad.top10, aes(x=reorder(fname, rndlasso), y=rndlasso, fill=short)) + geom_bar(stat = 'identity') +
  ylab('Importance') + 
  geom_hline(yintercept = 0.5, color='red') + coord_flip() + theme_Publication() + scale_fill_Publication1() +
  theme(legend.title = element_blank(), axis.title.y = element_blank())
ggsave('./figures/ICGC_ESAD.fs.top10.png', width = 8, height = 3, dpi=800)

# ICGC_LICA
lica = read.table('./figures/data/ICGC_LICA.select.tsv', header=T, sep='\t')
lica = merge(lica, meta, by.x='fname', by.y='feature_name')
lica = lica[order(lica$rndlasso, decreasing = T),]
lica.top10 = lica[1:10,]
lica.top10[c(2,4), 'short'] = 'Liver'
lica.top10[c(3,6), 'short'] = 'HepG2'
lica.top10[9, 'short'] = 'A549'
ggplot(lica.top10, aes(x=reorder(fname, rndlasso), y=rndlasso, fill=short)) + geom_bar(stat = 'identity') +
  ylab('Importance') + 
  geom_hline(yintercept = 0.5, color='red') + coord_flip() + theme_Publication() + scale_fill_Publication1() +
  theme(legend.title = element_blank(), axis.title.y = element_blank())
ggsave('./figures/ICGC_LICA.fs.top10.png', width = 8, height = 3, dpi=800)

# ICGC_PACA
paca = read.table('./figures/data/ICGC_PACA.select.tsv', header=T, sep='\t')
paca = merge(paca, meta, by.x='fname', by.y='feature_name')
paca = paca[order(paca$rndlasso, decreasing = T),]
paca.top10 = paca[1:10,]
paca.top10[5,'short'] = 'conservation'
paca.top10[10,'short'] = 'HepG2'
ggplot(paca.top10, aes(x=reorder(fname, rndlasso), y=rndlasso, fill=short)) + geom_bar(stat = 'identity') +
  ylab('Importance') + 
  geom_hline(yintercept = 0.5, color='red') + coord_flip() + theme_Publication() + scale_fill_Publication1() +
  theme(legend.title = element_blank(), axis.title.y = element_blank())
ggsave('./figures/ICGC_PACA.fs.top10.png', width = 8, height = 3, dpi=800)

# ICGA_BRCA
brca = read.table('./figures/data/ICGC_BRCA.select.tsv', header=T, sep='\t')
brca = merge(brca, meta, by.x='fname', by.y='feature_name')
brca = brca[order(brca$rndlasso, decreasing = T),]
brca.top10 = brca[1:10,]
brca.top10[5,'short'] = 'conservation'
brca.top10[6,'short'] = 'MCF7'
brca.top10[7,'short'] = 'HepG2'
brca.top10[9,'short'] = 'HeLa-S3'
ggplot(brca.top10, aes(x=reorder(fname, rndlasso), y=rndlasso, fill=short)) + geom_bar(stat = 'identity') +
  ylab('Importance') + 
  geom_hline(yintercept = 0.5, color='red') + coord_flip() + theme_Publication() + scale_fill_Publication1() +
  theme(legend.title = element_blank(), axis.title.y = element_blank())
ggsave('./figures/ICGC_BRCA.fs.top10.png', width = 8, height = 3, dpi=800)
