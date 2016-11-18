qqunif.plot<-function(pvector, title="Quantile-quantile plot of p-values") {
  # Thanks to Daniel Shriner at NHGRI for providing this code for creating expected and observed values
  o = -log10(sort(pvector,decreasing=F))
  e = -log10( 1:length(o)/length(o) )
  # you could use base graphics
  #plot(e,o,pch=19,cex=0.25, xlab=expression(Expected~~-log[10](italic(p))), ylab=expression(Observed~~-log[10](italic(p))),xlim=c(0,max(e)),ylim=c(0,max(e)))
  #lines(e,e,col="red")
  #You'll need ggplot2 installed to do the rest
  plot=qplot(e,o, xlim=c(0,max(e)), ylim=c(0,max(o))) + geom_abline(intercept=0,slope=1, col="red")
  plot=plot+ggtitle(title)
  plot=plot+scale_x_continuous(name=expression(Expected~~-log[10](italic(p))))
  plot=plot+scale_y_continuous(name=expression(Observed~~-log[10](italic(p))))
  plot
}
