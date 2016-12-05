# available CDS results
setwd('~/Desktop/DriverPower/results/CDS')
# Activedriver2 (ad2)
ad2 = read.table('ActiveDriver2.cds.tsv', header=T)
# compositeDriver (cd)
cd.res = read.table('compositeDriver.cds.tsv', header=T)
# dNdScv_NBR
nbr = read.table('dNdScv_NBR.cds.NBR.tsv', header=T)
dnds = read.table('dNdScv_NBR.cds.dNdScv.tsv', header=T)
# ExInAtor
ea.res = read.table('ExInAtor.cds.tsv', header=T)
# mutsig
mutsig = read.table('MutSig.cds.tsv', header=T)
# ncdDetect
ncdd = read.table('ncdDetect.cds.tsv', header=T)
# ncDriver
ncd = read.table('ncDriver.cds.tsv', header=T)
# oncodriveFML
odf.cadd = read.table('oncodriveFML.cds.cadd.tsv', header=T)
odf.vest3 = read.table('oncodriveFML.cds.vest3.tsv', header=T)

# check tumor types


# choose dnds and vest3 for more hits
combined = rbind(ad2, cd.res, dnds, ea.res, mutsig, ncdd, ncd, odf.vest3)
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

write.table(combined, './combined.cds.no.driverpower.tsv', sep='\t', row.names = F, quote=F)
