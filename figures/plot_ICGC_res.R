setwd('~/Desktop/DriverPower/figures/')
library(ggplot2)
library(ggrepel)
library(data.table)
source('./ggplot_theme.R')
plot_one_res = function(res.path, out.path, trim=TRUE){
  tb = read.table(res.path, header=T, sep='\t')
  if(trim){
    tb = tb[tb$p.value!=1,]
  }
  print(nrow(tb))
  print(nrow(tb[tb$q.value<0.1,]))
  tb['o'] = -log10(tb$p.value)
  tb['e'] = -log10(1:nrow(tb)/nrow(tb))
  tb['col'] = 'black'
  tb[tb$q.value<=0.1, 'col'] = "#386cb0"
  tb['text'] = tstrsplit(tb$element_ID, '::')[[3]]
  tb[tb$o>30, 'o'] = 30
  p = ggplot(tb, aes(e, o)) + geom_point(size=1, color=tb$col) + 
    geom_abline(intercept=0,slope=1, col="red") +
    geom_text_repel(data=tb[tb$q.value<0.1,], aes(e, o, label=text), color="#386cb0") +
    xlab(expression(Expected~~-log[10](italic(p)))) +
    ylab(expression(Observed~~-log[10](italic(p)))) +
    theme_Publication()
  ggsave(out.path, p, height = 6, width = 6)
}

# ICGC_LICA
name = c('cds', 'promCore', '3utr', '5utr', 'lncrna.ncrna', 'lncrna.promCore', 'enhancers')
for(i in name){
  plot_one_res(paste0('./data/ICGC_res/ICGC_LICA.', i, '.DriverPower.independent.txt'),
               paste0('./ICGC_figures/ICGC_LICA.', i, '.qqplot.png'), FALSE)
}

# ICGC_PACA
name = c('cds', 'promCore', '3utr', '5utr', 'lncrna.ncrna', 'lncrna.promCore', 'enhancers')
for(i in name){
  plot_one_res(paste0('./data/ICGC_res/ICGC_PACA.', i, '.DriverPower.independent.txt'),
               paste0('./ICGC_figures/ICGC_PACA.', i, '.qqplot.png'), FALSE)
}
