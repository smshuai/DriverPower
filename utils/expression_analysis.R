library(data.table)
# RNA
rna = fread('~/Downloads/tophat_star_fpkm_uq.v2_aliquot_gl_ncg.tsv',
            header=TRUE, sep='\t')
rnameta = read.table('~/Downloads/rnaseq_metadata.lite.with.barcode.tsv',
                     header=TRUE, sep='\t', check.names = FALSE)
# CN
cn = read.table('~/Downloads/all_samples.consensus_CN.by_gene.170214.txt',
                header=TRUE, sep='\t', check.names = FALSE)
cnv = read.table('~/Downloads/all_samples.consensus_level_calls.by_gene.170214.txt',
                 header=TRUE, sep='\t', check.names = FALSE)
# Purity
purity = read.table('~/Downloads/consensus.20170217.purity.ploidy.txt.gz', header=TRUE, sep='\t')
rnameta = merge(rnameta, purity, by.x='tumour_barcode', by.y='samplename')


library(ggplot2)
library(ggsci)
source('~/Desktop/DriverPower/utils/ggplot_theme.R')
df_one_gene = function(donors.in, rna, rnameta, cnv, gene.name='WDR74', gene.id='ENSG00000133316.11'){
  # obtain CNV
  g.cnv = cnv[cnv$`Gene Symbol`==gene.name, 4:ncol(cnv)]
  if(nrow(g.cnv) != 1) {
    warning('CNV is not right with ', nrow(g.cnv), ' row(s)')
    g.cnv[1,] = 1
  }
  g.cnv = data.frame(t(g.cnv), colnames(g.cnv))
  colnames(g.cnv) = c('cnv', 'tumour_barcode')
  # obtain expression
  g.rna = rna[grepl(gene.id, rna$feature), -1]
  if(nrow(g.rna) != 1) stop('RNA is not right with ', nrow(g.rna), ' row(s)')
  g.rna = data.frame(t(g.rna), colnames(g.rna))
  colnames(g.rna) = c('FPKM.UQ', 'aliquot_id')
  dat = merge(g.rna, rnameta)
  dat = merge(dat, g.cnv)
  dat['is_mut'] = ifelse(dat$icgc_donor_id %in% donors.in, 'MUT', 'WT')
  dat$is_mut[dat$is_tumour=='no'] = 'Normal'
  ylabel = "FPKM-UQ expression (log2 scale)"
  dat$is_mut = factor(dat$is_mut, levels = c('Normal', 'WT', 'MUT'))
  xlabs <- paste0(levels(dat$is_mut),"\n(N=",table(dat$is_mut),")")
  mod = glm(FPKM.UQ ~ is_mut + cnv + purity, data = dat[dat$is_tumour=='yes',], family = quasipoisson())
  res = anova(mod, test = 'LRT')
  cat('PanCan LRT pval =', res$`Pr(>Chi)`[2], '\n')
  p = qplot(x=is_mut, y=log2(FPKM.UQ+1), data=dat, geom = 'boxplot', fill=is_mut) + 
    ylab(ylabel) + ggtitle(paste0(gene.name, " - PanCan_No_Melanoma_Lymph")) + theme_Publication() + scale_fill_npg() +
    scale_x_discrete(labels=xlabs) + guides(fill=FALSE) + theme(axis.title.x = element_blank()) +
    annotate('text', x=-Inf, y=Inf, hjust=-.5, vjust=3,
             label=paste0('LRT p-value: ', as.character(round(res$`Pr(>Chi)`[2], 2))))
  print(p)
  # for tumours with >= 3 MUT
  ct = as.data.frame(table(dat$is_mut, dat$project))
  projects = as.character(ct[ct$Var1=='MUT' & ct$Freq>=3, 'Var2'])
  for (tumour in projects) {
    dat.p = subset(dat, project == tumour)
    mod = glm(FPKM.UQ ~ is_mut + cnv + purity, data = dat.p[dat.p$is_tumour=='yes',], family = quasipoisson())
    res = anova(mod, test = 'LRT')    
    cat(tumour, 'LRT pval = ', res$`Pr(>Chi)`[2], '\n')
    xlabs <- paste0(levels(dat.p$is_mut),"\n(N=",table(dat.p$is_mut),")")
    p = qplot(x=is_mut, y=log2(FPKM.UQ+1), data=dat.p, geom = 'boxplot', fill=is_mut) + 
      ylab(ylabel) + ggtitle(paste0(gene.name, " - ", tumour)) + theme_Publication() + scale_fill_npg() +
      scale_x_discrete(labels=xlabs, drop=FALSE) + guides(fill=FALSE) +
      theme(axis.title.x = element_blank()) +
      annotate('text', x=-Inf, y=Inf, hjust=-.5, vjust=3,
               label=paste0('LRT p-value: ', as.character(round(res$`Pr(>Chi)`[2], 2))))
    print(p)

  }
}

# example 1: WDR74
donors.in = scan('~/Downloads/donors.in.txt', what = character())
df_one_gene(donors.in, rna, rnameta, cnv, 'WDR74', 'ENSG00000133316.11')
df_one_gene(donors.in, rna, rnameta, cnv, 'ZFP36L2', 'ENSG00000152518.5')
df_one_gene(donors.in, rna, rnameta, cnv, 'RN7SK', 'ENSG00000202198.1')
df_one_gene(donors.in, rna, rnameta, cnv, 'FOXP2', 'ENSG00000128573.18')
df_one_gene(donors.in, rna, rnameta, cnv, 'GATA2', 'ENSG00000179348.7')

df_one_gene(donors.in, rna, rnameta, cnv, 'RNU1-2', 'ENSG00000207005.1')
df_one_gene(donors.in, rna, rnameta, cnv, 'CROCC', 'ENSG00000058453.12')

df_one_gene(donors.in, rna, rnameta, cnv, 'RNU6-9', 'ENSG00000207507.1')
df_one_gene(donors.in, rna, rnameta, cnv, 'MED16', 'ENSG00000175221.10')

df_one_gene(donors.in, rna, rnameta, cnv, 'CDK6', 'ENSG00000105810.5')

df_one_gene(donors.in, rna, rnameta, cnv, 'PARP16', 'ENSG00000138617.10')
