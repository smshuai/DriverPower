setwd('~/Desktop/DriverPower/figures/')
source('./ggplot_theme.R')
dat = read.table('./data/Liver-HCC.1mb.data.tsv',sep='\t', header = T)
X = dat[,4:ncol(dat)]
y = dat[,c(3, 2)] # N success, N all
y[,'fail'] = y$totcg - y$ct 
y = y[,-2]
# 500 test cases
test.idx = sample(1:2521, 500)
X.test = X[test.idx,]
X.train = X[-test.idx,]
y.test = as.matrix(y[test.idx,])
y.train = as.matrix(y[-test.idx,])

# glm
mod = glm(y.train ~ ., data = X.train, family = binomial())
mu.test = predict(mod, newdata = X.test, type = 'response')
mu.obs = y.test[,1] / rowSums(y.test)

mu.dat = data.frame(mu.test, mu.obs)
ggplot(mu.dat, aes(x=mu.obs, y=mu.test)) + geom_point() + xlab('Observed MR') +
  ylab('Predicted MR') + geom_abline(intercept=0,slope=1, color='red') +
  theme_Publication() + theme(axis.text.y = element_blank(), axis.text.x = element_blank())
ggsave('./pred.vs.obs.mu.1Mb.png', width = 2, height = 2, dpi=600)
