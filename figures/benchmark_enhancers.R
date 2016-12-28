setwd('~/DriverPower/results/')
source('../scripts/tune_func_cutoff.R')
source('../figures/ggplot_theme.R')
library(ggplot2)
library(ggrepel)

# singleTumor
# CADD
combined = read.table('./enhancers/sig.only/combined.enhancers.tsv', header=T, stringsAsFactors = F)
singleType = c('Liver-HCC', 'Panc-AdenoCA', 'Prost-AdenoCA', 'Breast-AdenoCA',
               'Kidney-RCC', 'CNS-Medullo', 'Ovary-AdenoCA', 'Skin-Melanoma', 'Lymph-BNHL', 'Eso-AdenoCA',
               'Lymph-CLL', 'CNS-PiloAstro', 'Panc-Endocrine', 'Stomach-AdenoCA', "Head-SCC",
               "ColoRect-AdenoCA", "Thy-AdenoCA", "Lung-SCC", "Uterus-AdenoCA",
               "Kidney-ChRCC", "Bone-Osteosarc", "CNS-GBM", 'Lung-AdenoCA', "Biliary-AdenoCA",
               'Bone-Leiomyo', 'Bladder-TCC', 'Myeloid-MPN', 'CNS-Oligo', 'Cervix-SCC')
file_names = grep('cadd', list.files('./enhancers/driverpower/', full.names=TRUE), value=TRUE)
keep = tolower(sapply(file_names, function(f) strsplit(basename(f), split = '.', fixed = TRUE)[[1]][1])) %in% tolower(singleType)
file_names = file_names[keep]
info = file.info(file_names) # remove empty files
file_names = rownames(info[info$size > 0, ])
func_cutoff = 70
dp = generate_dp(file_names, func_cutoff, ntest = 30816)
# fix name
dp$tumor = as.character(dp$tumor)
# fix Ca vs. CA
dp$tumor[dp$tumor == 'Breast-AdenoCa'] = 'Breast-AdenoCA'
dp$tumor[dp$tumor == 'Breast-LobularCa'] = 'Breast-LobularCA'
dp$tumor[dp$tumor == 'Eso-AdenoCa'] = 'Eso-AdenoCA'
dp$tumor[dp$tumor == 'Breast-AdenoCa'] = 'Breast-AdenoCA'
write.table(dp, './enhancers//sig.only/driverpower.enhancers.cadd70.singleTumor.tsv', sep = '\t', quote = F, row.names = F)
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
ggsave('../figures/plot/num.sig.singleType.enhancers.heatmap.png', heatmap.single, height = 10, width = 8)

combined.single = combined.single[combined.single$method != 'LARVA', ]
combined.single = combined.single[combined.single$tumor != 'Skin-Melanoma', ]
# combined.single = combined.single[combined.single$tumor != 'Lymph-BNHL', ]
dat.single = as.data.frame(table(as.data.frame(table(combined.single$id, combined.single$tumor))$Freq))
dat.single = dat.single[-1,]
bar.single = ggplot(dat.single, aes(x=Var1, y=Freq)) + geom_bar(stat = 'identity') + theme_Publication() + xlab('Number of Support')
ggsave('../figures/plot/nsup.singleType.enhancers.barplot.png', bar.single, height = 6, width = 6)

res = internal_benchmark(combined.single, cutoff.sup = 3, cutoff.cv = 2,
                         fig.path = '')

# EIGEN 93
file_names = grep('eigen', list.files('./enhancers/driverpower/', full.names=TRUE), value=TRUE)
keep = tolower(sapply(file_names, function(f) strsplit(basename(f), split = '.', fixed = TRUE)[[1]][1])) %in% tolower(singleType)
file_names = file_names[keep]
info = file.info(file_names) # remove empty files
file_names = rownames(info[info$size > 0, ])
func_cutoff = 75
dp = generate_dp(file_names, func_cutoff, ntest = 30816)
# fix name
dp$tumor = as.character(dp$tumor)
# fix Ca vs. CA
dp$tumor[dp$tumor == 'Breast-AdenoCa'] = 'Breast-AdenoCA'
dp$tumor[dp$tumor == 'Breast-LobularCa'] = 'Breast-LobularCA'
dp$tumor[dp$tumor == 'Eso-AdenoCa'] = 'Eso-AdenoCA'
write.table(dp, './enhancers/sig.only/driverpower.enhancers.eigen75.singleTumor.tsv', sep = '\t', quote = F, row.names = F)
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
ggsave('../figures/plot/num.sig.singleType.enhancers.heatmap.png', heatmap.single, height = 10, width = 8)

combined.single = combined.single[combined.single$method != 'LARVA', ]
combined.single = combined.single[combined.single$tumor != 'Skin-Melanoma', ]
dat.single = as.data.frame(table(as.data.frame(table(combined.single$id, combined.single$tumor))$Freq))
dat.single = dat.single[-1,]
bar.single = ggplot(dat.single, aes(x=Var1, y=Freq)) + geom_bar(stat = 'identity') + theme_Publication() + xlab('Number of Support')
ggsave('../figures/plot/nsup.singleType.enhancers.barplot.png', bar.single, height = 6, width = 6)

res = internal_benchmark(combined.single, cutoff.sup = 3, cutoff.cv = 2,
                         fig.path = '')
sub.res = combined.single
nsup = as.data.frame(table(sub.res$id, sub.res$tumor))
colnames(nsup) = c('id', 'tumor', 'num_sup')
sub.res = merge(sub.res, nsup)
sub.res$num_sup = ordered(sub.res$num_sup)
dp = sub.res[sub.res$method == 'DriverPower', ]
dp = dp[order(dp$tumor, dp$p), c(1,2,6)]
nsup$id = tstrsplit(nsup$id, '::')[[3]]
nsup[nsup$num_sup>1, ]
