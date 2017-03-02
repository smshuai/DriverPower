require(limma)
require(edgeR)
library(data.table)
library(ggplot2)

setwd('~/Desktop/DriverPower/figures/')
source('./ggplot_theme.R')
rnameta = read.table('./ALB/rnaseq_metadata.lite.tsv', header=T, sep='\t', stringsAsFactors = F)
rnameta = rnameta[rnameta$project_code %in% c('BLCA-US'),]
table(rnameta$is_tumour)

rnaseq = fread('~/Downloads/tophat2_count_extend.tsv', select = c('Feature', rnameta$aliquot_id))
rnaseq = as.data.frame(rnaseq)
row.names(rnaseq) = rnaseq$Feature
rnaseq.tumor = rnaseq[, c(rnameta$aliquot_id[rnameta$is_tumour=='yes'])]
condition <- colnames(count.dat)
donors.in = c('DO44097', "DO472", "DO477", "DO496",
              "DO498", "DO522", "DO555", "DO561",
              "DO669", 'DO695', 'DO804', 'DO822',
              'DO828', 'DO856')
condition = condition %in% rnameta$aliquot_id[rnameta$icgc_donor_id %in% donors.in]
# filter genes based on CPM
keep <- rowSums(rnaseq.tumor) > 23 
count.dat = rnaseq.tumor[keep,]
# create edgeR object
y <- DGEList(counts = count.dat, group = condition)
# normalization
y <- calcNormFactors(y)
y <- estimateDisp(y) # Disperison
et <- exactTest(y) # test
re.edgeR <- topTags(et, n=nrow(et))$table
# voom
v <- voom(y, plot = T)
design <- model.matrix(~ condition, data = as.data.frame(condition))
fit <- lmFit(v, design)
fit <- eBayes(fit)
summary(decideTests(fit))
res = as.matrix(topTable(fit, coef = 'conditionTRUE',  sort = "p"))
