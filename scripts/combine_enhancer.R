# available enhancer results
setwd('~/Desktop/DriverPower/results/enhancers/sig.only/')
# Activedriver2 (ad2)
ad2 = read.table('ActiveDriver2.enhancers.tsv', header=T, stringsAsFactors = F)
# dNdScv_NBR
nbr = read.table('NBR.enhancers.tsv', header=T, stringsAsFactors = F)
# mutsig
mutsig = read.table('./MutSig.enhancers.tsv', header=T, stringsAsFactors = F)
# ncdDetect
ncdd = read.table('./ncdDetect.enhancers.tsv', header=T, stringsAsFactors = F)
# regDriver
rd = read.table('./regDriver.enhancers.tsv', header=T, stringsAsFactors = F)
# oncodriveFML
odf = read.table('./oncodriveFML.enhancers.tsv', header=T, stringsAsFactors = F)
# LARVA
larva = read.table('./LARVA.enhancers.tsv', header=T, stringsAsFactors = F)
# ncDriver
ncd = read.table('./ncDriver.enhancers.tsv', header=T, stringsAsFactors = F)

combined = rbind(ad2, nbr, mutsig, ncdd, rd, odf, larva, ncd)

# unique tumor types before fix
sort(unique(combined$tumor)) # 71

# fix Ca vs. CA
combined$tumor[combined$tumor == 'Breast-AdenoCa'] = 'Breast-AdenoCA'
combined$tumor[combined$tumor == 'Breast-LobularCa'] = 'Breast-LobularCA'
combined$tumor[combined$tumor == 'Eso-AdenoCa'] = 'Eso-AdenoCA'
combined$tumor[combined$tumor == 'Breast-AdenoCa'] = 'Breast-AdenoCA'
combined$tumor[combined$tumor == 'Breast-AdenoCa'] = 'Breast-AdenoCA'


# fix meta names
combined$tumor[combined$tumor %in% c('Adenocarcinoma', "Adenocarcinoma-tumors", "meta_cohort_Adenocarcinoma_tumors")] = 'Adenocarcinoma_tumors'
combined$tumor[combined$tumor %in% c('Breast', 'Breast-tumors', "meta_cohort_Breast_tumors")] = 'Breast_tumors'
combined$tumor[combined$tumor %in% c('Carcinoma', "Carcinoma-tumors", "meta_cohort_Carcinoma_tumors")] = 'Carcinoma_tumors'
combined$tumor[combined$tumor %in% c('CNS', "CNS-tumors", "meta_cohort_CNS_tumors")] = 'CNS_tumors'
combined$tumor[combined$tumor %in% c('Digestive-tract', "Digestive-tract-tumors", 'Digestive_tract', "meta_cohort_Digestive_tract_tumors")] = 'Digestive_tract_tumors'
combined$tumor[combined$tumor %in% c('Female-reproductive-tract',
                                     'Female_reproductive_tract',
                                     "Female-reproductive-system-tumors",
                                     "meta_cohort_Female_reproductive_system_tumors")] = 'Female_reproductive_system_tumors'
combined$tumor[combined$tumor %in% c('Glioma', 'Glioma_tumors', "meta_cohort_Glioma_tumors")] = 'Glioma_tumors'
combined$tumor[combined$tumor %in% c('Hematopoietic-system', 'Hematopoietic_system',
                                     'Hematopoietic-system-tumors', "meta_cohort_Hematopoietic_tumors")] = 'Hematopoietic_tumors'
combined$tumor[combined$tumor %in% c('Kidney-tumors', 'Kidney', "meta_cohort_Kidney_tumors")] = 'Kidney_tumors'
combined$tumor[combined$tumor %in% c('Lung', 'Lung-tumors', "meta_cohort_Lung_tumors")] = 'Lung_tumors'
combined$tumor[combined$tumor %in% c('Lymphatic_system', 'Lymphatic-system', 'Lymph-tumors', "meta_cohort_Lymph_tumors")] = 'Lymph_tumors'
combined$tumor[combined$tumor %in% c('Sarcoma', 'Sarcoma-tumors', "meta_cohort_Sarcoma_tumors")] = 'Sarcoma_tumors'
combined$tumor[combined$tumor %in% c('Squamous', 'Squamous-tumors', "meta_cohort_Squamous_tumors")] = 'Squamous_tumors'
combined$tumor[combined$tumor %in% c('Myeloid', 'Myeloid-tumors', "meta_cohort_Myeloid_tumors")] = 'Myeloid_tumors'


# fix pancan names
combined$tumor[combined$tumor %in% c('All_cancers', 'pan', 'PANCANCER', 'PANCAN')] = 'PanCan'
combined$tumor[combined$tumor %in% c('All_cancers_no_Skin-Melanoma',
                                     'All_cancers-no-skin-melanoma',
                                     'PANCANCER_no_melanoma',
                                     'pan_noSkin-Melanoma')] = 'PanCan_No_Skin-Melanoma'
combined$tumor[combined$tumor %in% c('All_cancers-no-lymph', 'pan_noLymphatic_system')] = 'PanCan_No_Lymph_tumors'
combined$tumor[combined$tumor %in% c('All_cancers-no-skin-melanoma-lymph',
                                     'PANCANCER_no_melanoma_lymph', "PANCAN_noLymph_noSkin")] = 'PanCan_No_Melanoma_Lymph_tumors'


write.table(combined, './combined.enhancers.tsv', sep='\t', row.names = F, quote=F)

# create a heatmap
source('../../../figures/ggplot_theme.R')
library(ggplot2)
dat = as.data.frame(table(combined$tumor, combined$method))
ggplot(dat, aes(x=Var2, y=Var1, fill=Freq, label=Freq)) + geom_tile() + geom_text(color='white') + theme_Publication() +
  theme(legend.title = element_blank(), 
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(fill=F)

# meta by origin
metaOrigin = c('Adenocarcinoma_tumors', 'Glioma_tumors', 'Hematopoietic_tumors',
               'Sarcoma_tumors', 'Squamous_tumors', 'Carcinoma_tumors') # by origin
combined.origin = combined[combined$tumor %in% metaOrigin, ]
dat.origin = as.data.frame(table(combined.origin$tumor, combined.origin$method))
ggplot(dat.origin, aes(x=Var2, y=Var1, fill=Freq, label=Freq)) + geom_tile() + geom_text(color='white') + theme_Publication() +
  theme(legend.title = element_blank(), 
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(fill=F)
ggsave('../../../figures/plot/num.sig.metaOrigin.enhancers.heatmap.png')
combined.origin = combined.origin[combined.origin$tumor != 'Hematopoietic_tumors', ]
combined.origin = combined.origin[combined.origin$tumor != 'Carcinoma_tumors', ]
dat.origin = as.data.frame(table(as.data.frame(table(combined.origin$id, combined.origin$tumor))$Freq)[2:4])
ggplot(dat.origin, aes(x=Var1, y=Freq)) + geom_bar(stat = 'identity') + theme_Publication() + xlab('Number of Support')
ggsave('../../../figures/plot/nsup.enhancers.metaOrigin.barplot.png')

# meta by organ
metaOrgan = c('CNS_tumors', 'Digestive_tract_tumors', 'Female_reproductive_system_tumors',
              'Glioma_tumors', 'Kidney_tumors', 'Lung_tumors', 'Lymph_tumors', "Myeloid_tumors") # by organ
combined.organ = combined[combined$tumor %in% metaOrgan, ]
dat.organ = as.data.frame(table(combined.organ$tumor, combined.organ$method))
ggplot(dat.organ, aes(x=Var2, y=Var1, fill=Freq, label=Freq)) + geom_tile() + geom_text(color='white') + theme_Publication() +
  theme(legend.title = element_blank(), 
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(fill=F)
ggsave('../../../figures/plot/num.sig.metaOrgan.cds.heatmap.png')
combined.organ = combined.organ[combined.organ$tumor != 'Lymph_tumors', ]
dat.organ = as.data.frame(table(as.data.frame(table(combined.organ$id, combined.organ$tumor))$Freq)[2:8])
ggplot(dat.organ, aes(x=Var1, y=Freq)) + geom_bar(stat = 'identity') + theme_Publication() + xlab('Number of Support')
ggsave('../../../figures/plot/nsup.cds.metaOrgan.barplot.png')

# single Type
singleType = c(	'Liver-HCC', 'Panc-AdenoCA', 'Prost-AdenoCA', 'Breast-AdenoCA',
                'Kidney-RCC', 'CNS-Medullo', 'Ovary-AdenoCA', 'Skin-Melanoma', 'Lymph-BNHL', 'Eso-AdenoCA',
                'Lymph-CLL', 'CNS-PiloAstro', 'Panc-Endocrine', 'Stomach-AdenoCA', "Head-SCC",
                "ColoRect-AdenoCA", "Thy-AdenoCA", "Lung-SCC", "Uterus-AdenoCA",
                "Kidney-ChRCC", "Bone-Osteosarc", "CNS-GBM", 'Lung-AdenoCA', "Biliary-AdenoCA",
                'Bone-Leiomyo', 'Bladder-TCC', 'Myeloid-MPN', 'CNS-Oligo', 'Cervix-SCC')
combined.single = combined[combined$tumor %in% singleType, ]
dat.single = as.data.frame(table(combined.single$tumor, combined.single$method))
ggplot(dat.single, aes(x=Var2, y=Var1, fill=Freq, label=Freq)) + geom_tile() + geom_text(color='white') + theme_Publication() +
  theme(legend.title = element_blank(), 
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(fill=F)
ggsave('../../../figures/plot/num.sig.singleType.enhancers.heatmap.png')
combined.single = combined.single[combined.single$method != 'LARVA', ]
combined.single = combined.single[combined.single$method != 'ncDriver', ]
combined.single = combined.single[combined.single$tumor != 'Skin-Melanoma', ]
combined.single = combined.single[combined.single$tumor != 'Lymph-BNHL', ]

dat.single = as.data.frame(table(as.data.frame(table(combined.single$id, combined.single$tumor))$Freq)[2:4])
ggplot(dat.single, aes(x=Var1, y=Freq)) + geom_bar(stat = 'identity') + theme_Publication() + xlab('Number of Support')
ggsave('../../../figures/plot/nsup.enhancers.singleType.barplot.png')
