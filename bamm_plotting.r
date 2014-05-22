# Code for plotting results of BAMM analyses on phylogenies
# from 4 simulation scenarios
#
# 1) Pure niche conservatism, sim 3465
# 2) Energy gradient, sim 4065 (@ time = 100k)
# 3) Speciation gradient, sim 5525
# 4) Disturbance gradient, sim 3865 (need to run through 100k)

source('code/supplemental_analysis_functions.r')
source('../bamm-simulations/bamm_plotting.r')
require(ape)
require(BAMMtools)
require(RColorBrewer)


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

sims.to.load = c(4065, 3865, 5525)
sim.paths = c(#'Z:/SENCoutput/Hurlbert_and_Stegen_2014/raw_sim_output',
  'Z:/SENCoutput/Hurlbert_and_Stegen_2014/raw_sim_output',
  'Z:/SENCoutput/Hurlbert_and_Stegen_2014/disturb_sim_output',
  'Z:/Manuscripts/FrontiersTropicalDiversity/raw_sim_output')
burnInFrac = 0.2

if (0) { #this for loop not quite working yet
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


sim4065 = output.unzip('Z:/SENCoutput/Hurlbert_and_Stegen_2014/raw_sim_output', 4065)
sim3865 = output.unzip('Z:/SENCoutput/Hurlbert_and_Stegen_2014/disturb_sim_output', 3865)
sim5525 = output.unzip('Z:/Manuscripts/FrontiersTropicalDiversity/raw_sim_output', 5525)

all.pops4065 = sim4065$all.populations
all.pops3865 = sim3865$all.populations
all.pops5525 = sim5525$all.populations

extant.pops4065 = subset(all.pops4065, extant == 1)
extant.pops3865 = subset(all.pops3865, extant == 1)
extant.pops5525 = subset(all.pops5525, extant == 1)

phy4065 = sim4065$phylo.out
phy3865 = sim3865$phylo.out
phy5525 = sim5525$phylo.out

extant.phy4065 = read.tree('Z:/git/bamm-simulations/sim4065/extant_phy4065.tre')
extant.phy3865 = read.tree('Z:/git/bamm-simulations/sim3865/extant_phy3865.tre')
extant.phy5525 = read.tree('Z:/git/bamm-simulations/sim5525/extant_phy5525.tre')

edata4065 = getEventData(extant.phy4065, eventdata = "Z:/git/bamm-simulations/sim4065/sim4065_event_data.txt", burnin = burnInFrac)
edata3865 = getEventData(extant.phy3865, eventdata = "Z:/git/bamm-simulations/sim3865/sim3865_event_data.txt", burnin = burnInFrac)
edata5525 = getEventData(extant.phy5525, eventdata = "Z:/git/bamm-simulations/sim5525/sim5525_event_data.txt", burnin = burnInFrac)







