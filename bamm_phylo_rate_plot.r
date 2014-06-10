rates4065.30k = dtRates(edata4065.30k.3, 0.01)
x = rates4065.30k
colorbreaks <- assignColorBreaks2(x$dtrates$rates, 64, 's')

raterange = seq(min(x$dtrates$rates[[1]]), max(x$dtrates$rates[[1]]), length.out = 65)
colorobj = colorMap(raterange, pal, colorbreaks)


pdf('z:/manuscripts/frontierstropicaldiversity/bamm/rate_phylo_plot_4065-30k.pdf', height = 12, width = 6)
par(fig=c(0,1,0,0.2), oma = c(2,1,1,1))
plotRateThroughTime(edata4065.30k.3, axis.labels=F)
points(rep(31000,65), seq(min(x$dtrates$rates[[1]]), max(x$dtrates$rates[[1]]), length.out = 65), pch = 15, col = colorobj$cols)

par(fig=c(0,1,0.1,1), new=T)
bamm.plot(extant.phy4065.30k, edata4065.30k.3, extant.pops4065.30k)
mtext("Time before present", 1, outer=T)
dev.off()
