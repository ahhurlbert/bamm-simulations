# Code for plotting results of BAMM analyses on phylogenies
# from 4 simulation scenarios
#
# 1) Pure niche conservatism, sim 3465
# 2) Energy gradient, sim 4065 (@ time = 100k)
# 3) Speciation gradient, sim 5525
# 4) Disturbance gradient, sim 3865 (need to run through 100k)

source('../species-energy-simulation/code/supplemental_analysis_functions.r')
require(ape)
require(BAMMtools)
require(RColorBrewer)


bamm.plot = function(tree, edata, extantPops, colorbreaks = NULL, offset = 0.03) {
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
  
  xx = plot.bammdata(edata, lwd=1, colorbreaks = colorbreaks)  
  cols = colorRampPalette(c('violet', 'white','darkgreen'))(10)
  tiplabels(pch = 15, adj = (0.5 + offset), cex = 2, col = cols[home.regions$region])
  return(xx)
}

if (0) { #this for loop not quite working yet
sims.to.load = c(4065, 3865, 5525)
sim.paths = c(#'Z:/SENCoutput/Hurlbert_and_Stegen_2014/raw_sim_output',
  'Z:/SENCoutput/Hurlbert_and_Stegen_2014/raw_sim_output',
  'Z:/SENCoutput/Hurlbert_and_Stegen_2014/disturb_sim_output',
  'Z:/Manuscripts/FrontiersTropicalDiversity/raw_sim_output')
  for (i in 1:4) {
    simID = sims.to.load[i]
    tempsim = output.unzip(sim.paths[i], simID)
    assign(paste("sim", simID, sep = ""), tempsim)
    assign(paste("all.pops", simID, sep = ""), tempsim$all.populations)
    extant.pops = tempsim$all.populations[tempsim$all.populations$extant == 1,]
    assign(paste("extant.pops", simID, sep = ""), extant.pops)
    assign(paste("phy", simID, sep = ""), tempsim$phylo.out)
    extant.phy = drop.tip(tempsim$phylo.out, 
                          tempsim$phylo.out$tip.label[!tempsim$phylo.out$tip.label %in% extant.pops$spp.name])
    assign(paste("extant_phy", simID, sep = ""), extant.phy)
    assign(paste("edata", simID, sep = ""),
           getEventData(extant.phy, burnin = burnInFrac, eventdata = 
                          paste("Z:/git/bamm-simulations/sim", simID, "/sim", simID, "_event_data.txt", sep = "")))
    rm(list = c('tempsim', 'extant.pops', 'extant.phy'))
  }
} # end if


sim4065 = output.unzip('..species-energy-simulation/raw_sim_output', 4065)
sim3865 = output.unzip('..species-energy-simulation/raw_sim_output', 3865)
sim5525 = output.unzip('..species-energy-simulation/raw_sim_output', 5525)
sim3465 = output.unzip('..species-energy-simulation/raw_sim_output', 3465)

all.pops4065 = sim4065$all.populations
all.pops3865 = sim3865$all.populations
all.pops5525 = sim5525$all.populations
all.pops3465 = sim3465$all.populations

extant.pops4065 = subset(all.pops4065, extant == 1)
extant.pops3865 = subset(all.pops3865, extant == 1)
extant.pops5525 = subset(all.pops5525, extant == 1)
extant.pops3465 = subset(all.pops3465, extant == 1)
extant.pops4065.30k = subset(all.pops4065, time.of.origin <= 30000 & time.of.extinction > 30000)
extant.pops5525.30k = subset(all.pops5525, time.of.origin <= 30000 & time.of.extinction > 30000)

extant.phy4065 = read.tree('sim4065/extant_phy4065.tre')
extant.phy3865 = read.tree('sim3865/extant_phy3865.tre')
extant.phy5525 = read.tree('sim5525/extant_phy5525.tre')
extant.phy3465 = read.tree('sim3465/extant_phy3465.tre')
extant.phy4065.30k = read.tree('sim4065-30k/extant_phy4065_30k.tre')
extant.phy5525.30k = read.tree('sim5525-30k/extant_phy5525_30k.tre')

burnInFrac = 0.2
edata4065 = getEventData(extant.phy4065, eventdata = "sim4065/sim4065_event_data.txt", burnin = burnInFrac)
edata3865 = getEventData(extant.phy3865, eventdata = "sim3865/sim3865_event_data.txt", burnin = burnInFrac)
edata5525 = getEventData(extant.phy5525, eventdata = "sim5525/sim5525_event_data.txt", burnin = burnInFrac)
edata3465 = getEventData(extant.phy3465, eventdata = "sim3465/sim3465_event_data.txt", burnin = 0.05)
edata4065.30k = getEventData(extant.phy4065.30k, eventdata = "sim4065-30k/sim4065_30k_event_data.txt", burnin = burnInFrac)
edata5525.30k = getEventData(extant.phy5525.30k, eventdata = "sim5525-30k/sim5525_30k_event_data.txt", burnin = burnInFrac)


#Plotting some of these sims at t=30k and t=100k
#Need to make sure they are using the same color legend for comparison
pdf('bammplot_4_scenarios.pdf', height = 8, width = 10)
par(mfrow = c(1,2))
#Page 1
yy = bamm.plot(extant.phy5525.30k, edata5525.30k, extant.pops5525.30k)
mtext("Sim 5525, Speciation gradient, t = 30k", 3)
bamm.plot(extant.phy5525, edata5525, extant.pops5525, yy$colorbreaks)
mtext("Sim 5525, Speciation gradient, t = 100k", 3)
#Page 2
xx = bamm.plot(extant.phy4065.30k, edata4065.30k, extant.pops4065.30k)
mtext("Sim 4065, Energy gradient, t = 30k", 3)
bamm.plot(extant.phy4065, edata4065, extant.pops4065, xx$colorbreaks)
mtext("Sim 4065, Energy gradient, t = 100k", 3)
#Page 3
bamm.plot(extant.phy3865, edata3865, extant.pops3865)
mtext("Sim 3865, Disturbance gradient, t = 30k", 3)
#bamm.plot(extant.phy3465, edata3465, extant.pops3465)
#mtext("Sim 3465, Niche conservatism", 3)

dev.off()


