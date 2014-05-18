require(RColorBrewer)

burnInFrac = 0.2
tree5525 = read.tree('extant_phy5525.tre')
edata5525 <- getEventData(tree5525, eventdata = paste("event_data.txt", sep = ''), 
                      burnin = burnInFrac)

bamm.plot = function(simID, tree, edata, extantPops) {
  spp = unique(extantPops$spp.name)
  extantPops$env.diff = abs(extantPops$env.opt - extantPops$reg.env)
  home.reg = c()
  for (s in spp) {
    temp = subset(extantPops, spp.name == s)
    home.reg = rbind(home.reg, c(s, temp$region[temp$env.diff == max(temp$env.diff)]))
  }
  home.reg = as.data.frame(home.reg)
  names(home.reg) = c('spp.name', 'region')
  home.reg$spp.name = as.character(home.reg$spp.name)
  
  tips = data.frame(spp.name = as.character(tree$tip.label))
  home.regions = merge(tips, home.reg, by = 'spp.name', sort = F)
  
  plot.bammdata(edata, lwd=1)  
  cols = colorRampPalette(c('violet', 'white','darkgreen'))(10)
  tiplabels(pch = 15, adj = 0.52, cex = 2, col = cols[home.regions$region])
}

pdf(')
bamm.plot(4065, tree4065, edata4065, extant.pops4065)
mtext("Sim 4065, Energy Gradient, time = 100k")

