library(ggplot2)
setwd('~/DriverPower/results/CDS/')
# combined
combined = read.table('combined.cds.no.driverpower.tsv', header=T, stringsAsFactors = F)
# driverpower
dp = read.table('driverpower.cds.lasso.gmean.eigen.90.tsv', header=T, stringsAsFactors = F)
dp$method = 'DriverPower'
# tumor
# tumor = c("ColoRect-AdenoCA", "Eso-AdenoCa", "Liver-HCC", "Lung-AdenoCA", "Lung-SCC", "Panc-AdenoCA", "Stomach-AdenoCA")
# tumor = sort(unique(dp$tumor))
# tumor = tumor[tumor!='Carcinoma_tumors']
tumor = c('CNS_tumors', 'Digestive_tract_tumors', 'Female_reproductive_system_tumors',
          'Glioma_tumors', 'Kidney_tumors', 'Lung_tumors', 
          'Lymph_tumors') # by organ
tumor = c('Adenocarcinoma_tumors', 'Glioma_tumors', 'Hematopoietic_tumors', 'Sarcoma_tumors', 'Squamous_tumors') # by origin

sub.res = combined[combined$tumor %in% tumor, ]
sub.res = rbind(sub.res, dp[dp$tumor %in% tumor, ])




# number of support
nsup = as.data.frame(table(sub.res$id, sub.res$tumor))
colnames(nsup) = c('id', 'tumor', 'num_sup')
sub.res = merge(sub.res, nsup)
sub.res$num_sup = ordered(sub.res$num_sup)
# number of sig. regions per method
ggplot(sub.res, aes(x=method, fill=num_sup)) + geom_bar()
ggplot(sub.res, aes(x=method, fill=num_sup)) + geom_bar(position = 'fill')
# cross-validation
all_method = unique(sub.res$method)
cutoff = 2
res = data.frame(method = character(),
                 positive = numeric(),
                 true_positive = numeric(), 
                 false_positive = numeric(), 
                 precision = numeric(),
                 recall = numeric(),
                 F1 = numeric(), stringsAsFactors = FALSE)
for (i in 1:length(all_method)) {
  m = all_method[i]
  out.res = sub.res[sub.res$method != m, ]
  in.res  = sub.res[sub.res$method == m, ]
  ncv = as.data.frame(table(out.res$id, out.res$tumor))
  colnames(ncv) = c('id', 'tumor','num_cv')
  p = sum(ncv$num_cv >= cutoff)
  in.res = merge(in.res, ncv)
  tp = sum(in.res$num_cv >= cutoff)
  obs_false = sum(in.res$num_cv < cutoff)
  tpr = tp / p
  ppv = tp / nrow(in.res)
  f1 = 2*tpr*ppv/(ppv+tpr)
  res[i, 'method'] = m
  res[i, 2:ncol(res)] = c(p, tp, nrow(in.res)-tp, ppv, tpr, f1)
}
p = ggplot(res, aes(precision, recall, label=method)) + geom_point() + geom_text(hjust = 1, nudge_x = -0.01)
p + geom_text(data = res, aes(precision, recall, label=round(F1, 3)), hjust=0, nudge_x = 0.01) +
  xlim(0.55, 1.05)

