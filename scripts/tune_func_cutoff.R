library(ggplot2)


generate_dp = function(file_paths, cutoff){
  method = 'DriverPower'
  ans = data.frame(id = character(),
                   p = numeric(),
                   q = numeric(),
                   method = character(),
                   tumor = character(),
                   stringsAsFactors = FALSE)
  tot_sig = 0
  for (f in file_paths) {
    tumor = strsplit(basename(f), split = '.', fixed = TRUE)[[1]][1]
    cat(f, '\n')
    model.res = read.table(f, sep='\t', header=T, stringsAsFactors = F)
    one = cbind(tune_func_cutoff(cutoff, model.res), method, tumor)
    ans = rbind(ans, one)
    tot_sig = tot_sig + nrow(one)
  }
  cat('cutoff = ', cutoff, 'total sig. = ', tot_sig, '\n')
  return(ans)
}

tune_func_cutoff = function(cutoff, model.res){
  threshold = quantile(model.res$fscore, cutoff/100)
  
  MuAdj = model.res$Mu * threshold / model.res$fscore
  MuAdj[MuAdj > 1] = 1
  Pval = sapply(1:nrow(model.res),
                function(i) binom.test(model.res$nMutSample[i], model.res$Length[i], MuAdj[i], alternative = 'greater')$p.value)
  Qval = p.adjust(Pval, 'BH')
  # cat("cutoff =", cutoff, "Num_sig =", sum(Qval <= 0.1), "\n")
  ans = data.frame(id = model.res$binID[Qval <= 0.1],
                   p = Pval[Qval <= 0.1],
                   q = Qval[Qval <= 0.1], 
                   stringsAsFactors = FALSE)
  return(ans)
}

pr_table <- function(tumors, combined, file_names, func_cutoffs, positive_cutoff){
  ## generate precision recall table for different func_cutoff
  
  # extract data by tumor
  sub.res = combined[combined$tumor %in% tumors, ]
  # assume file name is tumor.XXXXX
  keep = sapply(file_names, function(f) strsplit(basename(f), split = '.', fixed = TRUE)[[1]][1]) %in% tumors
  file_paths = file_names[keep]
  # number of support
  nsup = as.data.frame(table(sub.res$id, sub.res$tumor))
  colnames(nsup) = c('id', 'tumor', 'num_sup')
  # Pseudo-Positive sets
  p = sum(nsup$num_sup >= positive_cutoff)
  # output format
  res = data.frame(cutoff = numeric(),
                   positive = numeric(),
                   true_positive = numeric(), 
                   false_positive = numeric(), 
                   precision = numeric(),
                   recall = numeric(),
                   F1 = numeric(), stringsAsFactors = FALSE)
  # main loop through cutoffs
  for (i in 1:length(func_cutoffs)) {
    cutoff = func_cutoffs[i]
    dp = generate_dp(file_paths, cutoff)
    dp = merge(dp, nsup, all.x = TRUE)
    dp$num_sup[is.na(dp$num_sup)] = 0
    tp = sum(dp$num_sup >= positive_cutoff)
    tpr = tp / p
    ppv = tp / nrow(dp)
    f1 = 2*tpr*ppv/(ppv+tpr)
    res[i, ] = c(cutoff, p, tp, nrow(dp)-tp, ppv, tpr, f1)
  }
  return(res)
}

internal_benchmark <- function(dp, combined, tumor, func_cutoff, cutoff.sup, cutoff.cv, fig.path){
  sub.res = combined[combined$tumor %in% tumor, ]
  sub.res = rbind(sub.res, dp[dp$tumor %in% tumor, ])
  # number of support
  nsup = as.data.frame(table(sub.res$id, sub.res$tumor))
  colnames(nsup) = c('id', 'tumor', 'num_sup')
  sub.res = merge(sub.res, nsup)
  sub.res$num_sup = ordered(sub.res$num_sup)
  # number of sig. regions per method
  barplot1 = ggplot(sub.res, aes(x=method, fill=num_sup)) + geom_bar() +
    theme_Publication() + scale_fill_Publication1() + 
    theme(legend.title = element_blank(), 
          axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1))
  barplot2 = ggplot(sub.res, aes(x=method, fill=num_sup)) + geom_bar(position = 'fill') +
    theme_Publication() + scale_fill_Publication1() +
    theme(legend.title = element_blank(), 
          axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1)) + ylab('Frequency')
  p.sup = sum(nsup$num_sup >= cutoff.sup)
  # cross-validation
  all_method = unique(sub.res$method)
  res = data.frame(method = character(),
                   positive_sup = numeric(),
                   true_positive_sup = numeric(), 
                   false_positive_sup = numeric(), 
                   precision_sup = numeric(),
                   recall_sup = numeric(),
                   F1_sup = numeric(),
                   positive_cv = numeric(),
                   true_positive_cv = numeric(), 
                   false_positive_cv = numeric(), 
                   precision_cv = numeric(),
                   recall_cv = numeric(),
                   F1_cv = numeric(),
                   stringsAsFactors = FALSE)
  for (i in 1:length(all_method)) {
    m = all_method[i]
    # res for other methods
    out.res = sub.res[sub.res$method != m, ]
    # res for benchmarked methods
    in.res  = sub.res[sub.res$method == m, ]
    ## benchmark with num_sup
    # True Positive
    tp.sup = sum(in.res$num_sup >= cutoff.sup)
    tpr.sup = tp.sup / p.sup
    ppv.sup = tp.sup / nrow(in.res)
    f1.sup = 2*tpr.sup*ppv.sup/(ppv.sup+tpr.sup)
    
    ## benchmark with num_cv
    ncv = as.data.frame(table(out.res$id, out.res$tumor))
    colnames(ncv) = c('id', 'tumor','num_cv')
    p.cv = sum(ncv$num_cv >= cutoff.cv)
    # keep all in.res rows
    in.res = merge(in.res, ncv, all.x = TRUE)
    in.res$num_cv[is.na(in.res$num_cv)] = 0
    # True Positive CV
    tp.cv = sum(in.res$num_cv >= cutoff.cv)
    # TPR (Recall) CV
    tpr.cv = tp.cv / p.cv
    # Precision CV
    ppv.cv = tp.cv / nrow(in.res)
    # F1 CV
    f1.cv = 2*tpr.cv*ppv.cv/(ppv.cv+tpr.cv)
    # output
    res[i, 'method'] = m
    res[i, 2:ncol(res)] = c(p.sup, tp.sup, nrow(in.res)-tp.sup, ppv.sup, tpr.sup, f1.sup,
                            p.cv, tp.cv, nrow(in.res)-tp.cv, ppv.cv, tpr.cv, f1.cv)
  }
  pr_plot.sup = ggplot(res, aes(recall_sup, precision_sup, color=method)) + geom_point() +
    geom_text_repel(label=paste(res$method, as.character(round(res$F1_sup, 3))), nudge_x = -0.005) +
    ylim(0.0, 1.0) + xlim(0.0, 1) +
    theme_Publication() + scale_color_Publication2() + guides(color=FALSE)
  pr_plot.cv = ggplot(res, aes(recall_cv, precision_cv, label=method, color=method)) + geom_point() +
    geom_text_repel(label=paste(res$method, as.character(round(res$F1_cv, 3))), nudge_x = -0.005) +
    ylim(0.0, 1.0) + xlim(0.0, 1) +
    theme_Publication() + scale_color_Publication2() + guides(color=FALSE)
  # combine plots
  png(fig.path, width = 960, height = 960)
  plot_list = list(barplot1, pr_plot.sup, barplot2, pr_plot.cv)
  multiplot(plotlist = plot_list, cols = 2)
  dev.off()
  return(res)
}

