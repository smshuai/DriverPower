library(ggplot2)
library(ggrepel)
setwd('~/Desktop/DriverPower/results/')
source('../scripts/tune_func_cutoff.R')
source('../figures/ggplot_theme.R')
##
## Generate data for precision-recall curve
##
# other methods result
combined = read.table('./CDS/sig.only/combined.cds.tsv', header=T, stringsAsFactors = F)
cutoffs = 80:99
# metaOrigin
metaOrigin = c('Adenocarcinoma_tumors', 'Glioma_tumors', 'Hematopoietic_tumors',
               'Sarcoma_tumors', 'Squamous_tumors') # by origin
# rndlasso + eigen + gmean
file_names = grep('cds.rndlasso.eigen.gmean', list.files('./CDS/driverpower/meta/', full.names=TRUE), value=TRUE)
resOrigin = pr_table(metaOrigin, combined, file_names, cutoffs, 3)
p.Origin = ggplot(resOrigin, aes(x=recall, y=precision, label=cutoff)) + geom_point() + geom_path(color="#386cb0") +
  geom_text_repel(nudge_x = 0.005, nudge_y = 0.005) + theme_Publication()
ggsave("../figures//plot/precision-recall.by.cutoff.metaOrigin.cds.rndlasso.eigen.gmean.png", p.Origin)
write.table(resOrigin, '../figures/data/precision-recall.metaOrigin.cds.rndlasso.eigen.gmean.3of5Method.tsv',
            sep='\t', row.names = F, quote = F)

# metaOrgan
metaOrgan = c('CNS_tumors', 'Digestive_tract_tumors', 'Female_reproductive_system_tumors',
              'Glioma_tumors', 'Kidney_tumors', 'Lung_tumors', 'Lymph_tumors') # by organ
# rndlasso + eigen + gmean
file_names = grep('cds.rndlasso.eigen.gmean', list.files('./CDS/driverpower/meta/', full.names=TRUE), value=TRUE)
resOrgan = pr_table(metaOrgan, combined, file_names, cutoffs, 3)
p.Organ = ggplot(resOrgan, aes(x=recall, y=precision, label=cutoff)) + geom_point() + geom_path(color="#386cb0") +
  geom_text_repel(nudge_x = 0.005, nudge_y = 0.005) + theme_Publication()
ggsave("../figures//plot/precision-recall.by.cutoff.metaOrgan.cds.rndlasso.eigen.gmean.png", p.Organ)
write.table(resOrgan, '../figures/data/precision-recall.metaOrgan.cds.rndlasso.eigen.gmean.3of5Method.tsv',
            sep='\t', row.names = F, quote = F)


##
## PromCore
##
# other methods result
combined = read.table('./promCore/sig.only/combined.promCore.tsv', header=T, stringsAsFactors = F)
cutoffs = 80:99
# metaOrigin
metaOrigin = c('Adenocarcinoma_tumors', 'Glioma_tumors', 'Hematopoietic_tumors',
               'Sarcoma_tumors', 'Squamous_tumors') # by origin
# rndlasso + eigen + gmean
file_names = grep('promCore.rndlasso.eigen.gmean', list.files('./promCore/driverpower/meta/', full.names=TRUE), value=TRUE)
resOrigin = pr_table(metaOrigin, combined, file_names, cutoffs, 3)
p.Origin = ggplot(resOrigin, aes(x=recall, y=precision, label=cutoff)) + geom_point() + geom_path(color="#386cb0") +
  geom_text_repel(nudge_x = 0.005, nudge_y = 0.005) + theme_Publication()
ggsave("../figures//plot/precision-recall.by.cutoff.metaOrigin.cds.rndlasso.eigen.gmean.png", p.Origin)
write.table(resOrigin, '../figures/data/precision-recall.metaOrigin.cds.rndlasso.eigen.gmean.3of5Method.tsv',
            sep='\t', row.names = F, quote = F)

# metaOrgan
metaOrgan = c('CNS_tumors', 'Digestive_tract_tumors', 'Female_reproductive_system_tumors',
              'Glioma_tumors', 'Kidney_tumors', 'Lung_tumors', 'Lymph_tumors') # by organ
# rndlasso + eigen + gmean
file_names = grep('promCore.rndlasso.eigen.gmean', list.files('./promCore/driverpower/meta/', full.names=TRUE), value=TRUE)
resOrgan = pr_table(metaOrgan, combined, file_names, cutoffs, 2)
p.Organ = ggplot(resOrgan, aes(x=recall, y=precision, label=cutoff)) + geom_point() + geom_path(color="#386cb0") +
  geom_text_repel(nudge_x = 0.005, nudge_y = 0.005) + theme_Publication()
ggsave("../figures//plot/precision-recall.by.cutoff.metaOrgan.cds.rndlasso.eigen.gmean.png", p.Organ)
write.table(resOrgan, '../figures/data/precision-recall.metaOrgan.cds.rndlasso.eigen.gmean.3of5Method.tsv',
            sep='\t', row.names = F, quote = F)