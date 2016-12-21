setwd('~/DriverPower/results/')
source('../scripts/tune_func_cutoff.R')
source('../figures/ggplot_theme.R')
library(ggplot2)
library(ggrepel)
##
## Internal Benchmark (promCore)
##
combined = read.table('./promCore/sig.only/combined.promCore.tsv', header=T, stringsAsFactors = F)
metaOrigin = c('Adenocarcinoma_tumors', 'Glioma_tumors', 'Hematopoietic_tumors',
               'Sarcoma_tumors', 'Squamous_tumors') # by origin
metaOrgan = c('CNS_tumors', 'Digestive_tract_tumors', 'Female_reproductive_system_tumors',
              'Glioma_tumors', 'Kidney_tumors', 'Lung_tumors', 'Lymph_tumors') # by organ

# metaOrgan
func_cutoff = 93
cutoff.sup = 3 # 4 out of 6
cutoff.cv = 2 # 3 out of 5
file_names = grep('promCore.rndlasso.eigen.gmean', list.files('./promCore/driverpower/meta/', full.names=TRUE), value=TRUE)
keep = sapply(file_names, function(f) strsplit(basename(f), split = '.', fixed = TRUE)[[1]][1]) %in% metaOrgan
file_paths = file_names[keep]
dp = generate_dp(file_paths, func_cutoff)
write.table(dp, './promCore/sig.only/driverpower.promCore.rndlasso.eigen93.gmean.metaOrgan.tsv', sep = '\t', quote = F, row.names = F)

res = internal_benchmark(dp, combined, metaOrgan, func_cutoff, cutoff.sup, cutoff.cv,
                         fig.path = '../figures/plot/benchmark.promCore.metaOrgan.func95.sup4.cv3.png')
write.table(res, '../figures/data/benchmark.cds.metaOrgan.func95.sup4.cv3.tsv', sep='\t', quote=F, row.names = F)
res = internal_benchmark(dp, combined, metaOrgan, func_cutoff, cutoff.sup=3, cutoff.cv=2,
                         fig.path = '../figures/plot/benchmark.cds.metaOrgan.func95.sup3.cv2.png')
write.table(res, '../figures/data/benchmark.cds.metaOrgan.func95.sup4.cv3.tsv', sep='\t', quote=F, row.names = F)

# metaOrigin
func_cutoff = 95
cutoff.sup = 4 # 4 out of 6
cutoff.cv = 3  # 3 out of 5
file_names = grep('cds.rndlasso.eigen.gmean', list.files('./CDS/driverpower/meta/', full.names=TRUE), value=TRUE)
keep = sapply(file_names, function(f) strsplit(basename(f), split = '.', fixed = TRUE)[[1]][1]) %in% metaOrigin
file_paths = file_names[keep]
dp = generate_dp(file_paths, func_cutoff)
res = internal_benchmark(dp, combined, metaOrigin, func_cutoff, cutoff.sup, cutoff.cv,
                         fig.path = '../figures/plot/benchmark.cds.metaOrigin.func95.sup4.cv3.png')
write.table(res, '../figures/data/benchmark.cds.metaOrigin.func95.sup4.cv3.tsv', sep='\t', quote=F, row.names = F)
res = internal_benchmark(dp, combined, metaOrigin, func_cutoff, cutoff.sup=3, cutoff.cv=2,
                         fig.path = '../figures/plot/benchmark.cds.metaOrigin.func95.sup3.cv2.png')
write.table(res, '../figures/data/benchmark.cds.metaOrigin.func95.sup3.cv2.tsv', sep='\t', quote=F, row.names = F)

# singleTumor
# CADD
singleType = c('Liver-HCC', 'Panc-AdenoCA', 'Prost-AdenoCA', 'Breast-AdenoCA',
               'Kidney-RCC', 'CNS-Medullo', 'Ovary-AdenoCA', 'Skin-Melanoma', 'Lymph-BNHL', 'Eso-AdenoCA',
               'Lymph-CLL', 'CNS-PiloAstro', 'Panc-Endocrine', 'Stomach-AdenoCA', "Head-SCC",
               "ColoRect-AdenoCA", "Thy-AdenoCA", "Lung-SCC", "Uterus-AdenoCA",
               "Kidney-ChRCC", "Bone-Osteosarc", "CNS-GBM", 'Lung-AdenoCA', "Biliary-AdenoCA",
               'Bone-Leiomyo', 'Bladder-TCC', 'Myeloid-MPN', 'CNS-Oligo', 'Cervix-SCC')
file_names = grep('cadd', list.files('./promCore/driverpower/', full.names=TRUE), value=TRUE)
keep = tolower(sapply(file_names, function(f) strsplit(basename(f), split = '.', fixed = TRUE)[[1]][1])) %in% tolower(singleType)
file_names = file_names[keep]
info = file.info(file_names) # remove empty files
file_names = rownames(info[info$size > 0, ])
func_cutoff = 40
dp = generate_dp(file_names, func_cutoff, ntest = 20164)
# fix name
dp$tumor = as.character(dp$tumor)
# fix Ca vs. CA
dp$tumor[dp$tumor == 'Breast-AdenoCa'] = 'Breast-AdenoCA'
dp$tumor[dp$tumor == 'Breast-LobularCa'] = 'Breast-LobularCA'
dp$tumor[dp$tumor == 'Eso-AdenoCa'] = 'Eso-AdenoCA'
dp$tumor[dp$tumor == 'Breast-AdenoCa'] = 'Breast-AdenoCA'
dp$tumor[dp$tumor == 'Breast-AdenoCa'] = 'Breast-AdenoCA'
write.table(dp, './promCore/sig.only/driverpower.promCore.cadd65.singleTumor.tsv', sep = '\t', quote = F, row.names = F)
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
ggsave('../figures/plot/num.sig.singleType.promCore.heatmap.png', heatmap.single, height = 10, width = 8)

combined.single = combined.single[combined.single$method != 'LARVA', ]
combined.single = combined.single[combined.single$method != 'ncDriver', ]
combined.single = combined.single[combined.single$tumor != 'Skin-Melanoma', ]
combined.single = combined.single[combined.single$tumor != 'Lymph-BNHL', ]
dat.single = as.data.frame(table(as.data.frame(table(combined.single$id, combined.single$tumor))$Freq))
dat.single = dat.single[-1,]
bar.single = ggplot(dat.single, aes(x=Var1, y=Freq)) + geom_bar(stat = 'identity') + theme_Publication() + xlab('Number of Support')
ggsave('../figures/plot/nsup.singleType.promCore.barplot.png', bar.single, height = 6, width = 6)

res = internal_benchmark(combined.single, cutoff.sup = 4, cutoff.cv = 3,
                         fig.path = '')

# EIGEN 93
file_names = grep('eigen', list.files('./promCore/driverpower/', full.names=TRUE), value=TRUE)
keep = tolower(sapply(file_names, function(f) strsplit(basename(f), split = '.', fixed = TRUE)[[1]][1])) %in% tolower(singleType)
file_names = file_names[keep]
info = file.info(file_names) # remove empty files
file_names = rownames(info[info$size > 0, ])
func_cutoff = 80
dp = generate_dp(file_names, func_cutoff, ntest = 20164)
# fix name
dp$tumor = as.character(dp$tumor)
# fix Ca vs. CA
dp$tumor[dp$tumor == 'Breast-AdenoCa'] = 'Breast-AdenoCA'
dp$tumor[dp$tumor == 'Breast-LobularCa'] = 'Breast-LobularCA'
dp$tumor[dp$tumor == 'Eso-AdenoCa'] = 'Eso-AdenoCA'
dp$tumor[dp$tumor == 'Breast-AdenoCa'] = 'Breast-AdenoCA'
dp$tumor[dp$tumor == 'Breast-AdenoCa'] = 'Breast-AdenoCA'
write.table(dp, './promCore/sig.only/driverpower.promCore.eigen80.singleTumor.tsv', sep = '\t', quote = F, row.names = F)
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
ggsave('../figures/plot/num.sig.singleType.promCore.heatmap.png', heatmap.single, height = 10, width = 8)

combined.single = combined.single[combined.single$method != 'LARVA', ]
combined.single = combined.single[combined.single$method != 'ncDriver', ]
combined.single = combined.single[combined.single$tumor != 'Skin-Melanoma', ]
combined.single = combined.single[combined.single$tumor != 'Lymph-BNHL', ]
dat.single = as.data.frame(table(as.data.frame(table(combined.single$id, combined.single$tumor))$Freq))
dat.single = dat.single[-1,]
bar.single = ggplot(dat.single, aes(x=Var1, y=Freq)) + geom_bar(stat = 'identity') + theme_Publication() + xlab('Number of Support')
ggsave('../figures/plot/nsup.singleType.promCore.barplot.png', bar.single, height = 6, width = 6)

res = internal_benchmark(combined.single, cutoff.sup = 4, cutoff.cv = 3,
                         fig.path = '')
sub.res = combined.single
nsup = as.data.frame(table(sub.res$id, sub.res$tumor))
colnames(nsup) = c('id', 'tumor', 'num_sup')
sub.res = merge(sub.res, nsup)
sub.res$num_sup = ordered(sub.res$num_sup)
dp = sub.res[sub.res$method == 'DriverPower', ]
dp$id = tstrsplit(dp$id, split = '::')[[3]]
dp = dp[order(dp$tumor, dp$p), c(1,2,6)]
