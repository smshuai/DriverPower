#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  args[2] = "plot.pdf"
} else if (length(args)>2) {
  stop('Unkown arguments supplied. Usage: plot_res.R res.txt plot.pdf')
}


# res = read.table('~/PanCan_No_Melanoma_Lymph.Cons100VertEl.DriverPower.res.tsv', header=TRUE, sep='\t')
res = read.table(args[1], header=TRUE, sep='\t')

# pval qqplot
source('~/PCAWG/DriverPower/utils/ggplot_theme.R')
pval_qqplot = function(pvals){
  require(ggplot2)
  # pvals can contain NA
  pvals = na.omit(pvals)
  observed <- sort(pvals)
  lobs <- -(log10(observed))
  expected <- c(1:length(observed)) 
  lexp <- -(log10(expected / (length(expected)+1)))
  p = qplot(x = lexp, y = lobs, geom = 'point') + 
    geom_abline(intercept = 0, slope = 1, color = 'red') +
    xlab(expression(Expected~~-log[10](italic(p)))) +
    ylab(expression(Observed~~-log[10](italic(p)))) +
    theme_Publication()
  return(p)
}

pdf(args[2])
pval_cols = colnames(res)[grepl('^p.', colnames(res))]
for (pcol in pval_cols) {
  p = pval_qqplot(res[,pcol])
  p = p + ggtitle(pcol)
  print(p)
}
dev.off()
