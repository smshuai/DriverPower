library(ggplot2)
library(ggthemes)
setwd('~/DriverPower/results/CDS/')
source('../scripts/tune_func_cutoff.R')

##
## Generate data for precision-recall curve
##
# other methods result
combined = read.table('./sig.only/combined.cds.no.driverpower.tsv', header=T, stringsAsFactors = F)
cutoffs = 80:99
# metaOrigin
metaOrigin = c('Adenocarcinoma_tumors', 'Glioma_tumors', 'Hematopoietic_tumors',
               'Sarcoma_tumors', 'Squamous_tumors') # by origin
# rndlasso + eigen + gmean
file_names = grep('cds.rndlasso.eigen.gmean', list.files('./driverpower/meta/', full.names=TRUE), value=TRUE)
resOrigin = pr_table(metaOrigin, combined, file_names, cutoffs, 3)
write.table(res, '../../figures/data/precision-recall.metaOrigin.cds.rndlasso.eigen.gmean.3of5Method.tsv',
            sep='\t', row.names = F, quote = F)
# metaOrgan
metaOrgan = c('CNS_tumors', 'Digestive_tract_tumors', 'Female_reproductive_system_tumors',
              'Glioma_tumors', 'Kidney_tumors', 'Lung_tumors', 'Lymph_tumors') # by organ
# rndlasso + eigen + gmean
file_names = grep('cds.rndlasso.eigen.gmean', list.files('./driverpower/meta/', full.names=TRUE), value=TRUE)
resOrgan = pr_table(metaOrgan, combined, file_names, cutoffs, 3)
write.table(res, '../../figures/data/precision-recall.metaOrgan.cds.rndlasso.eigen.gmean.3of5Method.tsv',
            sep='\t', row.names = F, quote = F)

# meta organ
## lasso cadd
organ.lc = read.delim('./precision-recall.metaOrgan.cds.lasso.cadd.gmean.3of5Method.tsv')
organ.lc = cbind(organ.lc, param='lasso+cadd')
## lasso eigen
organ.le = read.delim('./precision-recall.metaOrgan.cds.lasso.eigen.gmean.3of5Method.tsv')
organ.le = cbind(organ.le, param='lasso+eigen')
## rndlasso cadd
organ.rc = read.delim('./precision-recall.metaOrgan.cds.rndlasso.cadd.gmean.3of5Method.tsv')
organ.rc = cbind(organ.rc, param='rndlasso+cadd')
## rndlasso eigen
organ.re = read.delim('./precision-recall.metaOrgan.cds.rndlasso.eigen.gmean.3of5Method.tsv')
organ.re = cbind(organ.re, param='rndlasso+eigen')

dat = rbind(organ.lc, organ.le, organ.rc, organ.re)
ggplot(dat, aes(x=recall, y=precision, color=param)) + geom_line() + geom_path() + theme_Publication()


# metaOrgan rndlasso eigen
organ.re = read.delim('./precision-recall.metaOrgan.cds.rndlasso.eigen.gmean.3of5Method.tsv')
p.organ.re = ggplot(organ.re, aes(x=recall, y=precision, label=cutoff)) + geom_point() + geom_path(color="#386cb0") +
  geom_text(vjust=0, nudge_x = 0.005, nudge_y = 0.005) + theme_Publication()
ggsave("../plot/precision-recall.by.cutoff.metaOrgan.cds.rndlasso.eigen.gmean.png", p.organ.re)
# meta origin
origin.re = read.delim('./precision-recall.metaOrigin.cds.rndlasso.eigen.gmean.3of5Method.tsv')
p.origin.re = ggplot(origin.re, aes(x=recall, y=precision, label=cutoff)) + geom_point() + geom_path(color="#386cb0") +
  geom_text(vjust=0, nudge_x = 0.005, nudge_y = 0.005) + theme_Publication()
ggsave("../plot/precision-recall.by.cutoff.metaOrigin.cds.rndlasso.eigen.gmean.png", p.origin.re)

origin.le = read.delim('./precision-recall.metaOrigin.cds.lasso.eigen.gmean.3of5Method.tsv')
ggplot(origin.le, aes(x=recall, y=precision, label=cutoff)) + geom_point() + geom_path(color="#386cb0") +
  geom_text(vjust=0, nudge_x = 0.005, nudge_y = 0.005) + theme_Publication()
