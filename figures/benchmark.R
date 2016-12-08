setwd('~/DriverPower/results/')
source('../scripts/tune_func_cutoff.R')
source('../figures/ggplot_theme.R')
library(ggplot2)
library(ggrepel)


##
## Internal Benchmark (CDS)
##
combined = read.table('./CDS/sig.only/combined.cds.tsv', header=T, stringsAsFactors = F)
metaOrigin = c('Adenocarcinoma_tumors', 'Glioma_tumors', 'Hematopoietic_tumors',
               'Sarcoma_tumors', 'Squamous_tumors') # by origin
metaOrgan = c('CNS_tumors', 'Digestive_tract_tumors', 'Female_reproductive_system_tumors',
              'Glioma_tumors', 'Kidney_tumors', 'Lung_tumors', 'Lymph_tumors') # by organ

# metaOrgan
func_cutoff = 95
cutoff.sup = 4 # 4 out of 6
cutoff.cv = 3  # 3 out of 5
file_names = grep('cds.rndlasso.eigen.gmean', list.files('./CDS/driverpower/meta/', full.names=TRUE), value=TRUE)
keep = sapply(file_names, function(f) strsplit(basename(f), split = '.', fixed = TRUE)[[1]][1]) %in% metaOrgan
file_paths = file_names[keep]
dp = generate_dp(file_paths, func_cutoff)
res = internal_benchmark(dp, combined, metaOrgan, func_cutoff, cutoff.sup, cutoff.cv,
                         fig.path = '../figures/plot/benchmark.cds.metaOrgan.func95.sup4.cv3.png')
write.table(res, '../figures/data/benchmark.cds.metaOrgan.func95.sup4.cv3.tsv', sep='\t', quote=F, row.names = F)
res = internal_benchmark(dp, combined, metaOrgan, func_cutoff, cutoff.sup=3, cutoff.cv=2,
                         fig.path = '../figures/plot/benchmark.cds.metaOrgan.func95.sup3.cv2.png')
write.table(res, '../figures/data/benchmark.cds.metaOrgan.func95.sup4.cv3.tsv', sep='\t', quote=F, row.names = F)

# metaOrigin
func_cutoff = 95
cutoff.sup = 4 # 4 out of 6
cutoff.cv = 3  # 3 out of 5
file_names = grep('cds.rndlasso.eigen.gmean', list.files('./CDS/driverpower/meta/', full.names=TRUE), value=TRUE)
keep = sapply(file_names, function(f) strsplit(basename(f), split = '.', fixed = TRUE)[[1]][1]) %in% metaOrigin
file_paths = file_names[keep]
dp = generate_dp(file_paths, func_cutoff)
res = internal_benchmark(dp, combined, metaOrigin, func_cutoff, cutoff.sup, cutoff.cv,
                         fig.path = '../figures/plot/benchmark.cds.metaOrigin.func95.sup4.cv3.png')
write.table(res, '../figures/data/benchmark.cds.metaOrigin.func95.sup4.cv3.tsv', sep='\t', quote=F, row.names = F)
res = internal_benchmark(dp, combined, metaOrigin, func_cutoff, cutoff.sup=3, cutoff.cv=2,
                         fig.path = '../figures/plot/benchmark.cds.metaOrigin.func95.sup3.cv2.png')
write.table(res, '../figures/data/benchmark.cds.metaOrigin.func95.sup3.cv2.tsv', sep='\t', quote=F, row.names = F)
