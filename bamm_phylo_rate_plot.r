source('z:/git/species-energy-simulation/code/supplemental_analysis_functions.r')
library(BAMMtools)

assignColorBreaks2 = function (rates, NCOLORS = 64, spex = "s") 
{
  if (mode(rates) == "numeric") {
    bks <- quantile(rates, seq(0, 1, length.out = (NCOLORS + 
                                                     1)), na.rm = TRUE)
  }
  else if (mode(rates) == "list") {
    if (tolower(spex) == "s") {
      bks <- quantile(rates[[1]], seq(0, 1, length.out = (NCOLORS + 
                                                            1)), na.rm = TRUE)
    }
    else if (tolower(spex) == "e") {
      bks <- quantile(rates[[2]], seq(0, 1, length.out = (NCOLORS + 
                                                            1)), na.rm = TRUE)
    }
    else {
      bks <- quantile(rates[[1]] - rates[[2]], seq(0, 1, 
                                                   length.out = (NCOLORS + 1)), na.rm = TRUE)
    }
  }
  return(bks)
}


sim3465 = output.unzip('z:/sencoutput/hurlbert_and_stegen_2014/raw_sim_output/', 3465)
phy3465 = sim3465$phylo.out
all.pops3465 = sim3465$all.populations
max.time = max(all.pops3465$time.of.origin)
extant.pops3465 = subset(all.pops3465, extant==1)
extant.phy3465 = drop.tip(phy3465, phy3465$tip.label[!phy3465$tip.label %in% extant.pops3465$spp.name])
edata3465 = getEventData(extant.phy3465, 'sim3465/run1/sim3465_event_data.txt', burnin = 0.2)


rates4065.30k = dtRates(edata4065.30k.3, 0.01)
x = rates4065.30k
colorbreaks <- assignColorBreaks2(x$dtrates$rates, 64, 's')
pal = 'RdYlBu'

raterange = seq(min(x$dtrates$rates[[1]], na.rm = T), max(x$dtrates$rates[[1]], na.rm = T), length.out = 65)
colorobj = colorMap(raterange, pal, colorbreaks)


pdf('z:/manuscripts/frontierstropicaldiversity/bamm/rate_phylo_plot_3465k.pdf', height = 12, width = 6)
par(fig=c(0,1,0,0.25), oma = c(2,1,1,1))
plotRateThroughTime(edata3465, axis.labels = F, yticks = 4)
points(rep(1.03*max.time, 65), raterange, pch = 15, col = colorobj$cols)

par(fig=c(0,1,0.1,1), new=T)
bamm.plot(extant.phy4065.30k, edata4065.30k.3, extant.pops4065.30k)
mtext("Time before present", 1, outer=T)
dev.off()
