# available CDS results
setwd('~/DriverPower/results/promCore/sig.only/')
# Activedriver2 (ad2)
ad2 = read.table('ActiveDriver2.promCore.tsv', header=T, stringsAsFactors = F)
# compositeDriver (cd)
cd.res = read.table('compositeDriver.promCore.tsv', header=T, stringsAsFactors = F)
# dNdScv_NBR
nbr = read.table('dNdScv_NBR.promCore.tsv', header=T, stringsAsFactors = F)
# mutsig
mutsig = read.table('MutSig.promCore.tsv', header=T, stringsAsFactors = F)
# ncdDetect
ncdd = read.table('ncdDetect.promCore.tsv', header=T, stringsAsFactors = F)
# regDriver
rd = read.table('regDriver.promCore.tsv', header=T, stringsAsFactors = F)
# oncodriveFML
odf = read.table('oncodriveFML.promCore.tsv', header=T, stringsAsFactors = F)

combined = rbind(ad2, cd.res, nbr, mutsig, ncdd, rd, odf)
# fix Ca vs. CA
combined$tumor[combined$tumor == 'Breast-AdenoCa'] = 'Breast-AdenoCA'
combined$tumor[combined$tumor == 'Breast-LobularCa'] = 'Breast-LobularCA'
combined$tumor[combined$tumor == 'Eso-AdenoCa'] = 'Eso-AdenoCA'
combined$tumor[combined$tumor == 'Breast-AdenoCa'] = 'Breast-AdenoCA'
combined$tumor[combined$tumor == 'Breast-AdenoCa'] = 'Breast-AdenoCA'
# fix meta names
combined$tumor[combined$tumor %in% c('Adenocarcinoma', "Adenocarcinoma-tumors")] = 'Adenocarcinoma_tumors'
combined$tumor[combined$tumor %in% c('Breast', 'Breast-tumors')] = 'Breast_tumors'
combined$tumor[combined$tumor %in% c('Carcinoma', "Carcinoma-tumors")] = 'Carcinoma_tumors'
combined$tumor[combined$tumor %in% c('CNS', "CNS-tumors")] = 'CNS_tumors'
combined$tumor[combined$tumor %in% c('Digestive-tract', "Digestive-tract-tumors", 'Digestive_tract')] = 'Digestive_tract_tumors'
combined$tumor[combined$tumor %in% c('Female-reproductive-tract',
                                     'Female_reproductive_tract',
                                     "Female-reproductive-system-tumors")] = 'Female_reproductive_system_tumors'
combined$tumor[combined$tumor %in% c('Glioma', 'Glioma_tumors')] = 'Glioma_tumors'
combined$tumor[combined$tumor %in% c('Hematopoietic-system', 'Hematopoietic_system',
                                     'Hematopoietic-system-tumors')] = 'Hematopoietic_tumors'
combined$tumor[combined$tumor %in% c('Kidney-tumors', 'Kidney')] = 'Kidney_tumors'
combined$tumor[combined$tumor %in% c('Lung', 'Lung-tumors')] = 'Lung_tumors'
combined$tumor[combined$tumor %in% c('Lymphatic_system', 'Lymphatic-system', 'Lymph-tumors')] = 'Lymph_tumors'
combined$tumor[combined$tumor %in% c('Sarcoma', 'Sarcoma-tumors')] = 'Sarcoma_tumors'
combined$tumor[combined$tumor %in% c('Squamous', 'Squamous-tumors')] = 'Squamous_tumors'
combined$tumor[combined$tumor %in% c('Myeloid', 'Myeloid-tumors')] = 'Myeloid_tumors'
# fix pancan names
combined$tumor[combined$tumor %in% c('All_cancers', 'pan', 'PANCANCER', 'All-tumors')] = 'PanCan'
combined$tumor[combined$tumor %in% c('All_cancers_no_Skin-Melanoma',
                                     'All_cancers-no-skin-melanoma',
                                     'PANCANCER_no_melanoma',
                                     'pan_noSkin-Melanoma',
                                     'All-tumors-without-melanoma',
                                     "pancancer_no_melanoma")] = 'PanCan_No_Skin-Melanoma'
combined$tumor[combined$tumor %in% c('All_cancers-no-lymph',
                                     'pan_noLymphatic_system',
                                     "All-tumors-without-Lymphatic-system-tumors")] = 'PanCan_No_Lymph_tumors'
combined$tumor[combined$tumor %in% c('All_cancers-no-skin-melanoma-lymph',
                                     'PANCANCER_no_melanoma_lymph')] = 'PanCan_No_Melanoma_Lymph_tumors'

write.table(combined, './combined.promCore.tsv', sep='\t', row.names = F, quote=F)
