setwd('~/Desktop/DriverPower/results/CDS')


generate_dp = function(file_paths, cutoff){
  method = 'DriverPower'
  ans = data.frame(id = character(),
                   p = numeric(),
                   q = numeric(),
                   method = character(),
                   tumor = character(),
                   stringsAsFactors = FALSE)
  for (f in file_names) {
    tumor = strsplit(basename(f), split = '.', fixed = TRUE)[[1]][1]
    cat(f, '\n')
    model.res = read.table(f, sep='\t', header=T, stringsAsFactors = F)
    one = cbind(tune_func_cutoff(cutoff, model.res), method, tumor)
    ans = rbind(ans, one)
  }
  return(ans)
}

tune_func_cutoff = function(cutoff, model.res){
  threshold = quantile(model.res$fscore, cutoff/100)
  
  MuAdj = model.res$Mu * threshold / model.res$fscore
  MuAdj[MuAdj > 1] = 1
  Pval = sapply(1:nrow(model.res),
                function(i) binom.test(model.res$nMutSample[i], model.res$Length[i], MuAdj[i], alternative = 'greater')$p.value)
  Qval = p.adjust(Pval, 'BH')
  cat("cutoff =", cutoff, "Num_sig =", sum(Qval <= 0.1), "\n")
  ans = data.frame(id = model.res$binID[Qval <= 0.1],
                   p = Pval[Qval <= 0.1],
                   q = Qval[Qval <= 0.1], 
                   stringsAsFactors = FALSE)
  return(ans)
}

# other methods result
combined = read.table('~/final_consensus/analysis/cds/combined.cds.no.driverpower.tsv', header=T, stringsAsFactors = F)
tumor = c('Adenocarcinoma_tumors', 'Glioma_tumors', 'Hematopoietic_tumors', 'Sarcoma_tumors', 'Squamous_tumors') # by origin
# tumor = c('CNS_tumors', 'Digestive_tract_tumors', 'Female_reproductive_system_tumors',
#          'Glioma_tumors', 'Kidney_tumors', 'Lung_tumors', 
#          'Lymph_tumors') # by organ
sub.res = combined[combined$tumor %in% tumor, ]
# number of support
nsup = as.data.frame(table(sub.res$id, sub.res$tumor))
colnames(nsup) = c('id', 'tumor', 'num_sup')
# Pseudo-Positive sets
p = sum(nsup$num_sup >= 3)
cutoffs = 80:99
res = data.frame(cutoff = numeric(),
                 positive = numeric(),
                 true_positive = numeric(), 
                 false_positive = numeric(), 
                 precision = numeric(),
                 recall = numeric(),
                 F1 = numeric(), stringsAsFactors = FALSE)
# file_paths
file_names = grep('cds.rndlasso.eigen.gmean', list.files('./driverpower/meta/', full.names=TRUE), value=TRUE)
keep = sapply(file_names, function(f) strsplit(basename(f), split = '.', fixed = TRUE)[[1]][1]) %in% tumor
file_names = file_names[keep]
for (i in 1:length(cutoffs)) {
  cutoff = cutoffs[i]
  dp = generate_dp(file_names, cutoff)
  dp = merge(dp, nsup, all.x = TRUE)
  dp$num_sup[is.na(dp$num_sup)] = 0
  tp = sum(dp$num_sup >= 3)
  tpr = tp / p
  ppv = tp / nrow(dp)
  f1 = 2*tpr*ppv/(ppv+tpr)
  res[i, ] = c(cutoff, p, tp, nrow(dp)-tp, ppv, tpr, f1)
}
write.table(res, '~/final_consensus/analysis/cds/precision-recall.metaOrigin.cds.rndlasso.eigen.gmean.3of5Method.tsv',
            sep='\t', row.names = F, quote = F)
