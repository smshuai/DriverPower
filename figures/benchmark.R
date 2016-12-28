setwd('~/DriverPower/results/')
source('../scripts/tune_func_cutoff.R')
source('../figures/ggplot_theme.R')
library(ggplot2)
library(ggrepel)


##
## Internal Benchmark (CDS)
##
combined = read.table('./CDS/sig.only/combined.cds.tsv', header=T, stringsAsFactors = F)
metaOrigin = c('Adenocarcinoma_tumors', 'Glioma_tumors', 'Hematopoietic_tumors',
               'Sarcoma_tumors', 'Squamous_tumors') # by origin


# metaOrgan
# CADD 95
func_cutoff = 95
metaOrgan = c('CNS_tumors', 'Digestive_tract_tumors', 'Female_reproductive_system_tumors',
              'Glioma_tumors', 'Kidney_tumors', 'Lung_tumors', 'Lymph_tumors', "Myeloid_tumors") # by organ
file_names = grep('cadd', list.files('./CDS/driverpower/', full.names=TRUE), value=TRUE)
keep = tolower(sapply(file_names, function(f) strsplit(basename(f), split = '.', fixed = TRUE)[[1]][1])) %in% tolower(metaOrgan)
file_names = file_names[keep]
dp = generate_dp(file_names, func_cutoff, ntest = 20185)
# fix name
dp$tumor = as.character(dp$tumor)
write.table(dp, './CDS/sig.only/driverpower.cds.cadd95.metaOrgan.tsv', sep = '\t', quote = F, row.names = F)
# combined single
combined.organ = combined[combined$tumor %in% metaOrgan, ]
combined.organ = rbind(combined.organ, dp)
dat.organ = as.data.frame(table(combined.organ$tumor, combined.organ$method))
dat.organ$Freq[dat.organ$Freq == 0] = NA
heatmap.organ = ggplot(dat.organ, aes(x=Var2, y=Var1, fill=Freq, label=Freq)) + geom_tile() + geom_text(color='white') + theme_Publication() +
  theme(legend.title = element_blank(), 
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(fill=F)
ggsave('../figures/plot/num.sig.metaOrgan.cds.heatmap.png', heatmap.organ, height = 10, width = 8)
combined.organ = combined.organ[combined.organ$tumor != 'Lymph_tumors', ]
dat.organ = as.data.frame(table(as.data.frame(table(combined.organ$id, combined.organ$tumor))$Freq))[-1, ]
bar.organ = ggplot(dat.organ, aes(x=Var1, y=Freq)) + geom_bar(stat = 'identity') + theme_Publication() + xlab('Number of Support')
ggsave('../figures/plot/nsup.metaOrgan.cds.barplot.png', bar.organ, height = 6, width = 6)

res = internal_benchmark(combined.organ, cutoff.sup = 4, cutoff.cv = 3,
                         fig.path = '')

# EIGEN 93
func_cutoff=93
file_names = grep('eigen', list.files('./CDS/driverpower/', full.names=TRUE), value=TRUE)
keep = tolower(sapply(file_names, function(f) strsplit(basename(f), split = '.', fixed = TRUE)[[1]][1])) %in% tolower(metaOrgan)
file_names = file_names[keep]
dp = generate_dp(file_names, func_cutoff, ntest = 20185)
# fix name
dp$tumor = as.character(dp$tumor)
write.table(dp, './CDS/sig.only/driverpower.cds.eigen93.metaOrgan.tsv', sep = '\t', quote = F, row.names = F)
# combined single
combined.organ = combined[combined$tumor %in% metaOrgan, ]
combined.organ = rbind(combined.organ, dp)
dat.organ = as.data.frame(table(combined.organ$tumor, combined.organ$method))
dat.organ$Freq[dat.organ$Freq == 0] = NA
heatmap.organ = ggplot(dat.organ, aes(x=Var2, y=Var1, fill=Freq, label=Freq)) + geom_tile() + geom_text(color='white') + theme_Publication() +
  theme(legend.title = element_blank(), 
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(fill=F)
ggsave('../figures/plot/num.sig.metaOrgan.cds.heatmap.png', heatmap.organ, height = 10, width = 8)
combined.organ = combined.organ[combined.organ$tumor != 'Lymph_tumors', ]
dat.organ = as.data.frame(table(as.data.frame(table(combined.organ$id, combined.organ$tumor))$Freq))[-1, ]
bar.organ = ggplot(dat.organ, aes(x=Var1, y=Freq)) + geom_bar(stat = 'identity') + theme_Publication() + xlab('Number of Support')
ggsave('../figures/plot/nsup.metaOrgan.cds.barplot.png', bar.organ, height = 6, width = 6)

res = internal_benchmark(combined.organ, cutoff.sup = 4, cutoff.cv = 3,
                         fig.path = '')

# metaOrigin
# CADD
func_cutoff = 95
metaOrigin = c('Adenocarcinoma_tumors', 'Glioma_tumors', 'Hematopoietic_tumors',
               'Sarcoma_tumors', 'Squamous_tumors', 'Carcinoma_tumors') # by origin
file_names = grep('cadd', list.files('./CDS/driverpower/', full.names=TRUE), value=TRUE)
keep = tolower(sapply(file_names, function(f) strsplit(basename(f), split = '.', fixed = TRUE)[[1]][1])) %in% tolower(metaOrigin)
file_names = file_names[keep]
dp = generate_dp(file_names, func_cutoff, ntest = 20185)
# fix name
dp$tumor = as.character(dp$tumor)
write.table(dp, './CDS/sig.only/driverpower.cds.cadd95.metaOrigin.tsv', sep = '\t', quote = F, row.names = F)
combined.origin = combined[combined$tumor %in% metaOrigin, ]
combined.origin = rbind(combined.origin, dp)
dat.origin = as.data.frame(table(combined.origin$tumor, combined.origin$method))
dat.origin[dat.origin == 0] = NA
heatmap.origin = ggplot(dat.origin, aes(x=Var2, y=Var1, fill=Freq, label=Freq)) + geom_tile() + geom_text(color='white') + theme_Publication() +
  theme(legend.title = element_blank(), 
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(fill=F)
ggsave('../figures/plot/num.sig.metaOrigin.cds.heatmap.png', heatmap.origin, height = 10, width = 8)
combined.origin = combined.origin[combined.origin$tumor != 'Hematopoietic_tumors', ]
combined.origin = combined.origin[combined.origin$tumor != 'Carcinoma_tumors', ]
dat.origin = as.data.frame(table(as.data.frame(table(combined.origin$id, combined.origin$tumor))$Freq))[-1,]
bar.origin = ggplot(dat.origin, aes(x=Var1, y=Freq)) + geom_bar(stat = 'identity') + theme_Publication() + xlab('Number of Support')
ggsave('../figures/plot/nsup.metaOrigin.cds.barplot.png', bar.origin,  height = 6, width = 6)
res = internal_benchmark(combined.origin, cutoff.sup = 4, cutoff.cv = 3,
                         fig.path = '')
# Eigen
func_cutoff = 95
metaOrigin = c('Adenocarcinoma_tumors', 'Glioma_tumors', 'Hematopoietic_tumors',
               'Sarcoma_tumors', 'Squamous_tumors', 'Carcinoma_tumors') # by origin
file_names = grep('eigen', list.files('./CDS/driverpower/', full.names=TRUE), value=TRUE)
keep = tolower(sapply(file_names, function(f) strsplit(basename(f), split = '.', fixed = TRUE)[[1]][1])) %in% tolower(metaOrigin)
file_names = file_names[keep]
dp = generate_dp(file_names, func_cutoff, ntest = 20185)
# fix name
dp$tumor = as.character(dp$tumor)
write.table(dp, './CDS/sig.only/driverpower.cds.eigen95.metaOrigin.tsv', sep = '\t', quote = F, row.names = F)
combined.origin = combined[combined$tumor %in% metaOrigin, ]
combined.origin = rbind(combined.origin, dp)
dat.origin = as.data.frame(table(combined.origin$tumor, combined.origin$method))
dat.origin[dat.origin == 0] = NA
heatmap.origin = ggplot(dat.origin, aes(x=Var2, y=Var1, fill=Freq, label=Freq)) + geom_tile() + geom_text(color='white') + theme_Publication() +
  theme(legend.title = element_blank(), 
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(fill=F)
ggsave('../figures/plot/num.sig.metaOrigin.cds.heatmap.png', heatmap.origin, height = 10, width = 8)
combined.origin = combined.origin[combined.origin$tumor != 'Hematopoietic_tumors', ]
combined.origin = combined.origin[combined.origin$tumor != 'Carcinoma_tumors', ]
dat.origin = as.data.frame(table(as.data.frame(table(combined.origin$id, combined.origin$tumor))$Freq))[-1,]
bar.origin = ggplot(dat.origin, aes(x=Var1, y=Freq)) + geom_bar(stat = 'identity') + theme_Publication() + xlab('Number of Support')
ggsave('../figures/plot/nsup.metaOrigin.cds.barplot.png', bar.origin,  height = 6, width = 6)
res = internal_benchmark(combined.origin, cutoff.sup = 4, cutoff.cv = 3,
                         fig.path = '')

# singleTumor
func_cutoff = 95
singleType = c('Liver-HCC', 'Panc-AdenoCA', 'Prost-AdenoCA', 'Breast-AdenoCA',
               'Kidney-RCC', 'CNS-Medullo', 'Ovary-AdenoCA', 'Skin-Melanoma', 'Lymph-BNHL', 'Eso-AdenoCA',
               'Lymph-CLL', 'CNS-PiloAstro', 'Panc-Endocrine', 'Stomach-AdenoCA', "Head-SCC",
               "ColoRect-AdenoCA", "Thy-AdenoCA", "Lung-SCC", "Uterus-AdenoCA",
               "Kidney-ChRCC", "Bone-Osteosarc", "CNS-GBM", 'Lung-AdenoCA', "Biliary-AdenoCA",
               'Bone-Leiomyo', 'Bladder-TCC', 'Myeloid-MPN', 'CNS-Oligo', 'Cervix-SCC')
# CADD
file_names = grep('cadd', list.files('./CDS/driverpower/', full.names=TRUE), value=TRUE)
keep = tolower(sapply(file_names, function(f) strsplit(basename(f), split = '.', fixed = TRUE)[[1]][1])) %in% tolower(singleType)
file_names = file_names[keep]
info = file.info(file_names) # remove empty files
file_names = rownames(info[info$size > 0, ])
dp = generate_dp(file_names, func_cutoff, ntest = 20185)
# fix name
dp$tumor = as.character(dp$tumor)
# fix Ca vs. CA
dp$tumor[dp$tumor == 'Breast-AdenoCa'] = 'Breast-AdenoCA'
dp$tumor[dp$tumor == 'Breast-LobularCa'] = 'Breast-LobularCA'
dp$tumor[dp$tumor == 'Eso-AdenoCa'] = 'Eso-AdenoCA'
dp$tumor[dp$tumor == 'Breast-AdenoCa'] = 'Breast-AdenoCA'
dp$tumor[dp$tumor == 'Breast-AdenoCa'] = 'Breast-AdenoCA'
write.table(dp, './promCore/sig.only/driverpower.promCore.cadd95.singleTumor.tsv', sep = '\t', quote = F, row.names = F)
# combined single
combined.single = combined[combined$tumor %in% singleType, ]
combined.single = rbind(combined.single, dp)
dat.single = as.data.frame(table(combined.single$tumor, combined.single$method))
dat.single$Freq[dat.single$Freq == 0] = NA
heatmap.single = ggplot(dat.single, aes(x=Var2, y=Var1, fill=Freq, label=Freq)) + geom_tile() + geom_text(color='white') + theme_Publication() +
  theme(legend.title = element_blank(), 
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(fill=F)
ggsave('../figures/plot/num.sig.singleType.cds.heatmap.png', heatmap.single, height = 10, width = 8)

combined.single = combined.single[combined.single$method != 'LARVA', ]
combined.single = combined.single[combined.single$method != 'ncDriver', ]
combined.single = combined.single[combined.single$tumor != 'Skin-Melanoma', ]
combined.single = combined.single[combined.single$tumor != 'Lymph-BNHL', ]
dat.single = as.data.frame(table(as.data.frame(table(combined.single$id, combined.single$tumor))$Freq))
dat.single = dat.single[-1,]
bar.single = ggplot(dat.single, aes(x=Var1, y=Freq)) + geom_bar(stat = 'identity') + theme_Publication() + xlab('Number of Support')
ggsave('../figures/plot/nsup.singleType.cds.barplot.png', bar.single, height = 6, width = 6)

res = internal_benchmark(combined.single, singleType, func_cutoff, cutoff.sup = 4, cutoff.cv = 3,
                         fig.path = '../figures/plot/benchmark.cds.singleType.cadd95.png')

# EIGEN 93
file_names = grep('eigen', list.files('./CDS/driverpower/', full.names=TRUE), value=TRUE)
keep = tolower(sapply(file_names, function(f) strsplit(basename(f), split = '.', fixed = TRUE)[[1]][1])) %in% tolower(singleType)
file_names = file_names[keep]
dp = generate_dp(file_names, func_cutoff, ntest = 20185)
# fix name
dp$tumor = as.character(dp$tumor)
# fix Ca vs. CA
dp$tumor[dp$tumor == 'Breast-AdenoCa'] = 'Breast-AdenoCA'
dp$tumor[dp$tumor == 'Breast-LobularCa'] = 'Breast-LobularCA'
dp$tumor[dp$tumor == 'Eso-AdenoCa'] = 'Eso-AdenoCA'
dp$tumor[dp$tumor == 'Breast-AdenoCa'] = 'Breast-AdenoCA'
dp$tumor[dp$tumor == 'Breast-AdenoCa'] = 'Breast-AdenoCA'
write.table(dp, './CDS/sig.only/driverpower.eigen93.singleTumor.tsv', sep = '\t', quote = F, row.names = F)
# combined single
combined.single = combined[combined$tumor %in% singleType, ]
combined.single = rbind(combined.single, dp)
dat.single = as.data.frame(table(combined.single$tumor, combined.single$method))
dat.single$Freq[dat.single$Freq == 0] = NA
heatmap.single = ggplot(dat.single, aes(x=Var2, y=Var1, fill=Freq, label=Freq)) + geom_tile() + geom_text(color='white') + theme_Publication() +
  theme(legend.title = element_blank(), 
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(fill=F)
ggsave('../figures/plot/num.sig.singleType.cds.heatmap.png', heatmap.single, height = 10, width = 8)

combined.single = combined.single[combined.single$method != 'LARVA', ]
combined.single = combined.single[combined.single$method != 'ncDriver', ]
combined.single = combined.single[combined.single$tumor != 'Skin-Melanoma', ]
combined.single = combined.single[combined.single$tumor != 'Lymph-BNHL', ]
dat.single = as.data.frame(table(as.data.frame(table(combined.single$id, combined.single$tumor))$Freq))
dat.single = dat.single[-1,]
bar.single = ggplot(dat.single, aes(x=Var1, y=Freq)) + geom_bar(stat = 'identity') + theme_Publication() + xlab('Number of Support')
ggsave('../figures/plot/nsup.singleType.cds.barplot.png', bar.single, height = 6, width = 6)

res = internal_benchmark(combined.single, singleType, func_cutoff, cutoff.sup = 4, cutoff.cv = 3,
                         fig.path = '../figures/plot/benchmark.cds.singleType.eigen93.png')


