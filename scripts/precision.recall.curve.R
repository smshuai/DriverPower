library(ggplot2)
library(ggthemes)
setwd('~/final_consensus/analysis/cds/')
source('./ggplot_theme.R')
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
ggplot(dat, aes(x=recall, y=precision, color=param)) + geom_line() + geom_point() 
ggplot(organ.re, aes(x=recall, y=precision, label=cutoff)) + geom_point() + geom_line(color="#386cb0") +
  geom_text(vjust=0, nudge_x = 0.005, nudge_y = 0.005) + theme_Publication()

