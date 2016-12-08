# available CDS results
setwd('~/DriverPower/results/CDS/sig.only/')
# Activedriver2 (ad2)
ad2 = read.table('ActiveDriver2.cds.tsv', header=T, stringsAsFactors = F)
# compositeDriver (cd)
cd.res = read.table('compositeDriver.cds.tsv', header=T, stringsAsFactors = F)
# dNdScv_NBR
nbr = read.table('dNdScv_NBR.cds.NBR.tsv', header=T, stringsAsFactors = F)
dnds = read.table('dNdScv_NBR.cds.dNdScv.tsv', header=T, stringsAsFactors = F)
# ExInAtor
ea.res = read.table('ExInAtor.cds.tsv', header=T, stringsAsFactors = F)
# mutsig
mutsig = read.table('MutSig.cds.tsv', header=T, stringsAsFactors = F)
# ncdDetect
ncdd = read.table('ncdDetect.cds.tsv', header=T, stringsAsFactors = F)
# ncDriver
ncd = read.table('ncDriver.cds.tsv', header=T, stringsAsFactors = F)
# oncodriveFML
odf.cadd = read.table('oncodriveFML.cds.cadd.tsv', header=T, stringsAsFactors = F)
odf.vest3 = read.table('oncodriveFML.cds.vest3.tsv', header=T, stringsAsFactors = F)

# check tumor types


# choose dnds and vest3 for more hits
combined = rbind(ad2, cd.res, dnds, ea.res, mutsig, ncdd, ncd, odf.cadd)

# fix Ca vs. CA
combined$tumor[combined$tumor == 'Breast-AdenoCa'] = 'Breast-AdenoCA'
combined$tumor[combined$tumor == 'Breast-LobularCa'] = 'Breast-LobularCA'
combined$tumor[combined$tumor == 'Eso-AdenoCa'] = 'Eso-AdenoCA'
combined$tumor[combined$tumor == 'Breast-AdenoCa'] = 'Breast-AdenoCA'
combined$tumor[combined$tumor == 'Breast-AdenoCa'] = 'Breast-AdenoCA'


# fix names
combined$tumor[combined$tumor == 'Adenocarcinoma'] = 'Adenocarcinoma_tumors'
combined$tumor[combined$tumor == 'Breast'] = 'Breast_tumors'
combined$tumor[combined$tumor == 'Carcinoma'] = 'Carcinoma_tumors'
combined$tumor[combined$tumor == 'CNS'] = 'CNS_tumors'
combined$tumor[combined$tumor == 'Digestive-tract'] = 'Digestive_tract_tumors'
combined$tumor[combined$tumor == 'Digestive_tract'] = 'Digestive_tract_tumors'
combined$tumor[combined$tumor == 'Female-reproductive-tract'] = 'Female_reproductive_system_tumors'
combined$tumor[combined$tumor == 'Female_reproductive_tract'] = 'Female_reproductive_system_tumors'
combined$tumor[combined$tumor == 'Glioma'] = 'Glioma_tumors'
combined$tumor[combined$tumor == 'Hematopoietic-system'] = 'Hematopoietic_tumors'
combined$tumor[combined$tumor == 'Hematopoietic_system'] = 'Hematopoietic_tumors'
combined$tumor[combined$tumor == 'Kidney'] = 'Kidney_tumors'
combined$tumor[combined$tumor == 'Lung'] = 'Lung_tumors'
combined$tumor[combined$tumor == 'Lymphatic_system'] = 'Lymph_tumors'
combined$tumor[combined$tumor == 'Lymphatic-system'] = 'Lymph_tumors'
combined$tumor[combined$tumor == 'Sarcoma'] = 'Sarcoma_tumors'
combined$tumor[combined$tumor == 'Squamous'] = 'Squamous_tumors'
combined$tumor[combined$tumor == 'Myeloid'] = 'Myeloid_tumors'
# fix pancan names
combined$tumor[combined$tumor %in% c('All_cancers', 'pan', 'PANCANCER')] = 'PanCan'
combined$tumor[combined$tumor %in% c('All_cancers_no_Skin-Melanoma',
                                     'All_cancers-no-skin-melanoma',
                                     'PANCANCER_no_melanoma',
                                     'pan_noSkin-Melanoma')] = 'PanCan_No_Skin-Melanoma'
combined$tumor[combined$tumor %in% c('All_cancers-no-lymph', 'pan_noLymphatic_system')] = 'PanCan_No_Lymph_tumors'
combined$tumor[combined$tumor %in% c('All_cancers-no-skin-melanoma-lymph',
                                     'PANCANCER_no_melanoma_lymph')] = 'PanCan_No_Melanoma_Lymph_tumors'



write.table(combined, './combined.cds.tsv', sep='\t', row.names = F, quote=F)
