library(ggplot2)
library(ggthemes)
setwd('~/DriverPower/figures/data/')
source('../ggplot_theme.R')
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
ggsave("../plot/precision-recall.by.cutoff.metaOrgan.cds.rndlasso.eigen.gmean.png", p.origin.re)
# meta origin
origin.re = read.delim('./precision-recall.metaOrigin.cds.rndlasso.eigen.gmean.3of5Method.tsv')
p.origin.re = ggplot(origin.re, aes(x=recall, y=precision, label=cutoff)) + geom_point() + geom_path(color="#386cb0") +
  geom_text(vjust=0, nudge_x = 0.005, nudge_y = 0.005) + theme_Publication()
ggsave("../plot/precision-recall.by.cutoff.metaOrigin.cds.rndlasso.eigen.gmean.png", p.origin.re)

origin.le = read.delim('./precision-recall.metaOrigin.cds.lasso.eigen.gmean.3of5Method.tsv')
ggplot(origin.le, aes(x=recall, y=precision, label=cutoff)) + geom_point() + geom_path(color="#386cb0") +
  geom_text(vjust=0, nudge_x = 0.005, nudge_y = 0.005) + theme_Publication()
