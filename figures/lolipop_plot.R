library(ggplot2)
# ALB
mut = read.table('./ALB/Liver-HCC.ALB.mut.tsv', stringsAsFactors = F, sep='\t')
bed = read.table('./ALB/gc19.exon.bed', stringsAsFactors = F, sep='\t')
plot_lollipop(mut, bed, "gc19_pc.cds::gencode::ALB::ENSG00000163631.12")
bed.promCore = read.table('./ALB/gc19.promCore.bed', stringsAsFactors = F, sep='\t')
plot_lollipop(mut, bed.promCore, 'gc19_pc.promCore::gencode::ALB::ENSG00000163631.12')
bed.3utr = read.table('./ALB/gc19.3utr.bed', stringsAsFactors = F, sep='\t')
plot_lollipop(mut, bed.3utr, 'gc19_pc.3utr::gencode::ALB::ENSG00000163631.12')
bed.ss = read.table('./ALB/gc19_pc.ss.5col.bed', stringsAsFactors = F, sep='\t')
plot_lollipop(mut, bed.ss, 'gc19_pc.ss::gencode::ALB::ENSG00000163631.12')

# GNAS
mut = read.table('./GNAS/Panc-AdenoCA.GNAS.CDS.mut.tsv', stringsAsFactors = F, sep='\t')
plot_lollipop(mut, bed, 'gc19_pc.cds::gencode::GNAS::ENSG00000087460.19')

# GPR126
mut = read.table('./GPR126/Bladder-TCC.gpr126.mut.tsv', stringsAsFactors = F, sep='\t')
bed.enhancer = read.table('./ALB/enhancers.5col.bed', stringsAsFactors = F, sep='\t')
plot_lollipop(mut, bed.enhancer, 'enhancers::chr6:142705600-142706400::NA::NA')
mut = read.table('./GPR126/Breast-AdenoCa.gpr126.mut.tsv', stringsAsFactors = F, sep='\t')
plot_lollipop(mut, bed.enhancer, 'enhancers::chr6:142705600-142706400::NA::NA')


plot_lollipop <-function(mut, bed, element_ID){
  element = subset(bed, V4==element_ID)
  element = element[order(element$V2),]
  element[,'relative_start'] = 0
  element[,'relative_end'] = 0
  cum_len = 0
  for (i in 1:nrow(element)){
    new_len = cum_len+element[i, 3]-element[i, 2]
    element$relative_start[i] = cum_len
    element$relative_end[i] = new_len - 1
    cum_len = new_len
  }
  start.eles = element$V2
  mut.in.ele = subset(mut, V13==element_ID)[, c('V2','V4','V14')]
  mut.in.ele = merge(mut.in.ele, element, by.x='V14', by.y='V5')
  mut.in.ele['coords'] = mut.in.ele$V2.x - mut.in.ele$V2.y + mut.in.ele$relative_start
  dat = as.data.frame(table(mut.in.ele[,c('V4.x', 'coords')]))
  dat = dat[dat$Freq>0,]
  dat$coords = as.numeric(as.character(dat$coords))
  p = ggplot() + geom_vline(xintercept = element$relative_start[-1], linetype="dotted") + theme_Publication() +
    geom_segment(data = dat, aes(x=coords, xend=coords, y=0, yend=Freq)) + xlab('Relative coordinates') +
    geom_point(data = dat, aes(x=coords, y=Freq, color=V4.x, size=10)) + ylim(-1,NA) +
    annotate('rect', xmin=0, xmax=max(element$relative_end), ymin=-0.5, ymax=0, fill='darkgrey') +
    guides(size=FALSE) + theme(legend.title = element_blank()) + ylab('Number of mutations')
  return(p)
}


# GNAS
