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
donors.in = scan('~/Downloads/donors.in.WDR74', what = character())
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
df_one_gene(donors.in, rna, rnameta, cnv, 'RNU5A-1', 'ENSG00000199568.1')
df_one_gene(donors.in, rna, rnameta, cnv, 'PARP16', 'ENSG00000138617.10')

# res
res = read.table('~/Downloads/Driver Regions - No-PCAWG set merged results.tsv',
                 stringsAsFactors = F, header=T, sep='\t')
targets = lapply(strsplit(res$PCAWG_elements, ','),
       function(x) {
         if (x[1] != '.'){
         parts = tstrsplit(x, '::')
         unique(paste(parts[[2]], parts[[3]], parts[[4]], sep = '::'))
         }else{
           x
         }
       }
)
targets.append = sapply(strsplit(res$PCAWG_elements, ','),
                        function(x) {
                          if (x[1] != '.'){
                            parts = tstrsplit(x, '::')
                            uniq_id = unique(parts[[3]])
                            uniq_id = uniq_id[!grepl('NA', uniq_id)]
                            paste(uniq_id, collapse =';')
                          }else{
                            x
                          }
                        })
targets.append[104] = 'GPR126'
res$Target = targets.append
# no target regions 
no_target = subset(res, Target == '.')
no_regions = sapply(strsplit(no_target$Element_ID, ','), function(x) unique(tstrsplit(x, '::')[[1]]))
no_regions = do.call(c, no_regions)
table(no_regions)

# uniq targets
uniq_targets = unique(do.call(c, targets))
# remove microRNA prom
uniq_targets = uniq_targets[!grepl('^chr', uniq_targets)]
# Add GPR126
uniq_targets = c(uniq_targets, 'gencode::GPR126::ENSG00000112414.10')
# 101 uniq targets in total

# Compare normal expression vs tumour
head(rnameta)
normal_vs_tumour <- function(g.name, g.id, rnameta, rna, projects){
  pvals = numeric(1+length(projects))
  names(pvals) = c('PanCan_No_Melanoma_Lymph', projects)
  g.rna = rna[rna$feature == g.id, -1]
  if(nrow(g.rna) != 1) stop('RNA is not right with ', nrow(g.rna), ' row(s)')
  g.rna = data.frame(t(g.rna), colnames(g.rna))
  colnames(g.rna) = c('FPKM.UQ', 'aliquot_id')
  dat = merge(g.rna, rnameta)
  pval = LRT_test(dat)
  pvals[1] = pval
  p = plot_dat(dat, pval)
  # pancan
  p = p + ggtitle(paste0(g.name, " - ", 'PanCan_No_Melanoma_Lymph'))
  print(p)
  for (tumour in projects){
    dat.p = subset(dat, project == tumour)
    pval = LRT_test(dat.p)
    pvals[tumour] = pval
    p = plot_dat(dat.p, pval)
    p = p + ggtitle(paste0(g.name, " - ", tumour))
    print(p)
  }
  return(pvals)
}

rnameta$purity[rnameta$is_tumour=='no']=0
LRT_test <- function(dat, form='FPKM.UQ ~ is_tumour + purity', plot_var='is_tumour'){
  # LRT test
  mod = glm(form, data = dat, family = quasipoisson())
  pval = anova(mod, test = 'LRT')$`Pr(>Chi)`[2]
  return(pval)
}
plot_dat <- function(dat, pval, plot_var='is_tumour'){
  dat['logFPKM'] = log2(dat$FPKM.UQ+1)
  dat[plot_var] = as.factor(dat[,plot_var])
  xlabs <- paste0(levels(dat[,plot_var]),"\n(N=",table(dat[plot_var]),")")
  p = ggplot(dat, aes_string(x=plot_var, y='logFPKM')) + geom_boxplot(outlier.shape = NA) +
    geom_jitter(aes_string(color=plot_var)) +
    ylab("FPKM-UQ expression (log2 scale)") + theme_Publication() + scale_color_npg() +
    scale_x_discrete(labels=xlabs) + guides(color=FALSE) + theme(axis.title.x = element_blank()) +
    annotate('text', x=-Inf, y=Inf, hjust=-.5, vjust=3,
             label=paste0('LRT p-value: ', as.character(round(pval, 4))))
  return(p)
}

ct = as.data.frame(table(rnameta$project, rnameta$is_tumour))
projects = as.character(ct$Var1[ct$Var2=='no' & ct$Freq>=2])
uniq_targets = uniq_targets[-25] # remove .
write.table(uniq_targets,'~/Downloads/uniq_targets.txt', sep='\t', quote = F, row.names = F)
uniq_targets = scan('~/Downloads/uniq_targets.txt', what = character())
pdf('~/Downloads/exp_normal_vs_tumour.pdf')
all_pvals = list()
for (gene in uniq_targets){
  gname = strsplit(gene, '::')[[1]][2]
  gid = paste0(strsplit(gene, '::')[[1]][1], '::', gname, '::tc_', strsplit(gene, '::')[[1]][3])
  print(gid)
  all_pvals[gid] = normal_vs_tumour(gname, gid, rnameta, rna, projects)
}
dev.off()
normal_vs_tumour('MIR663A', 'ENSG00000227195.4', rnameta, rna, projects)
