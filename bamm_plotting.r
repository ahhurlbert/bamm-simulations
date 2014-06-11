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

# Temporary functions replacing assignColorBreaks and plot.bammdata
# to remedy a problem with passing rates of NA to quantiles() within
# assignColorBreaks. Here, I've simply added a na.rm = TRUE argument.
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

# And now plot.bammdata2 simply calls assignColorBreaks2 instead of assignColorBreaks
# I've also added a lab.col argument to specify tip label color
# (I make it white to create marginal white space, what a hack)
plot.bammdata2 = function (x, method = "phylogram", vtheta = 5, rbf = 0.001, show = TRUE, 
                           labels = FALSE, legend = FALSE, spex = "s", lwd = 1, cex = 1, 
                           pal = "RdYlBu", colorbreaks = NULL, par.reset = TRUE, lab.col = 'black', ...) 
{
  if ("bammdata" %in% class(x)) {
    if (attributes(x)$order != "cladewise") {
      stop("Function requires tree in 'cladewise' order")
    }
    phy <- as.phylo.bammdata(x)
  }
  else stop("Object ephy must be of class bammdata")
  if (!is.binary.tree(phy)) {
    stop("Function requires fully bifurcating tree")
  }
  if (any(phy$edge.length == 0)) {
    warning("Tree contains zero length branches. Rates for these will be NA and coerced to zero")
  }
  if (!("dtrates" %in% names(x))) {
    x <- dtRates(x, 0.01)
  }
  if (is.null(colorbreaks)) {
    colorbreaks <- assignColorBreaks2(x$dtrates$rates, 64, 
                                      spex)
  }
  if (x$type == "trait") {
    if (sum(is.na(x$dtrates$rates))) {
      warning(sprintf("Found %d NA phenotypic rates. Coercing to zero.", 
                      sum(is.na(x$dtrates$rates))))
      x$dtrates$rates[is.na(x$dtrates$rates)] <- 0
    }
    colorobj <- colorMap(x$dtrates$rates, pal, colorbreaks)
  }
  else if (x$type == "diversification") {
    if (sum(is.na(x$dtrates$rates[[1]]))) {
      warning(sprintf("Found %d NA speciation rates. Coercing to zero.", 
                      sum(is.na(x$dtrates$rates[[1]]))))
      x$dtrates$rates[[1]][is.na(x$dtrates$rates[[1]])] <- 0
    }
    if (sum(is.na(x$dtrates$rates[[2]]))) {
      warning(sprintf("Found %d NA extinction rates. Coercing to zero.", 
                      sum(is.na(x$dtrates$rates[[2]]))))
      x$dtrates$rates[[2]][is.na(x$dtrates$rates[[2]])] <- 0
    }
    if (tolower(spex) == "s") {
      colorobj <- colorMap(x$dtrates$rates[[1]], pal, colorbreaks)
    }
    else if (tolower(spex) == "e") {
      colorobj <- colorMap(x$dtrates$rates[[2]], pal, colorbreaks)
    }
    else {
      colorobj <- colorMap(x$dtrates$rates[[1]] - x$dtrates$rates[[2]], 
                           pal, colorbreaks)
    }
  }
  else {
    stop("Unrecognized/corrupt bammdata class. Type does not equal 'trait' or 'diversification'")
  }
  edge.color <- colorobj$cols
  tH <- max(branching.times(phy))
  phy$begin <- x$begin
  phy$end <- x$end
  tau <- x$dtrates$tau
  if (method == "polar") {
    ret <- setPolarTreeCoords(phy, vtheta, rbf)
    rb <- tH * rbf
    p <- mkdtsegsPolar(ret$segs[-1, ], tau, x$edge)
  }
  else if (method == "phylogram") {
    ret <- setPhyloTreeCoords(phy)
    p <- mkdtsegsPhylo(ret$segs[-1, ], tau, x$edge)
  }
  else {
    stop("Unimplemented method")
  }
  x0 <- c(ret$segs[1, 1], p[, 1])
  x1 <- c(ret$segs[1, 3], p[, 2])
  y0 <- c(ret$segs[1, 2], p[, 3])
  y1 <- c(ret$segs[1, 4], p[, 4])
  offset <- table(p[, 5])[as.character(unique(p[, 5]))]
  arc.color <- c(edge.color[1], edge.color[match(unique(p[, 
                                                          5]), p[, 5]) + offset])
  edge.color <- c(edge.color[1], edge.color)
  if (show) {
    if (length(list(...))) {
      op <- par(no.readonly = TRUE)
      par(...)
    }
    plot.new()
    ofs <- 0
    if (labels) {
      ofs <- max(nchar(phy$tip.label) * 0.03 * cex)
    }
    if (method == "polar") {
      plot.window(xlim = c(-1, 1) + c(-rb, rb) + c(-ofs, 
                                                   ofs), ylim = c(-1, 1) + c(-rb, rb) + c(-ofs, 
                                                                                          ofs), asp = 1)
      segments(x0, y0, x1, y1, col = edge.color, lwd = lwd, 
               lend = 2)
      arc(0, 0, ret$arcs[, 1], ret$arcs[, 2], c(rb, rb + 
                                                  phy$end/tH), border = arc.color, lwd = lwd)
      if (labels) {
        for (k in 1:length(phy$tip.label)) {
          text(ret$segs[-1, ][phy$edge[, 2] == k, 3], 
               ret$segs[-1, ][phy$edge[, 2] == k, 4], phy$tip.label[k], 
               cex = cex, srt = (180/pi) * ret$arcs[-1, 
                                                    ][phy$edge[, 2] == k, 1], adj = c(0, NA),
               col = lab.col)
        }
      }
    }
    if (method == "phylogram") {
      plot.window(xlim = c(0, 1 + ofs), ylim = c(0, phy$Nnode * 
                                                   1/(phy$Nnode + 1)))
      segments(x0[-1], y0[-1], x1[-1], y1[-1], col = edge.color[-1], 
               lwd = lwd, lend = 2)
      isTip <- phy$edge[, 2] <= phy$Nnode + 1
      isTip <- c(FALSE, isTip)
      segments(ret$arcs[!isTip, 1], ret$arcs[!isTip, 2], 
               ret$arcs[!isTip, 3], ret$arcs[!isTip, 4], col = arc.color[!isTip], 
               lwd = lwd, lend = 2)
      if (labels) {
        text(ret$segs[-1, ][phy$edge[, 2] <= phy$Nnode + 
                              1, 3], ret$segs[-1, ][phy$edge[, 2] <= phy$Nnode + 
                                                      1, 4], phy$tip.label, cex = cex, pos = 4, offset = 0.25,
             col = lab.col)
      }
    }
    if (legend) {
      rateLegend(colorobj$colsdensity)
    }
  }
  index <- order(as.numeric(rownames(ret$segs)))
  if (method == "phylogram") {
    assign("last_plot.phylo", list(type = "phylogram", direction = "rightwards", 
                                   Ntip = phy$Nnode + 1, Nnode = phy$Nnode, edge = phy$edge, 
                                   xx = ret$segs[index, 3], yy = ret$segs[index, 4], 
                                   pp = par(no.readonly = TRUE)), envir = .PlotPhyloEnv)
  }
  else if (method == "polar") {
    assign("last_plot.phylo", list(type = "fan", Ntip = phy$Nnode + 
                                     1, Nnode = phy$Nnode, edge = phy$edge, xx = ret$segs[index, 
                                                                                          3], yy = ret$segs[index, 4], theta = ret$segs[index, 
                                                                                                                                        5], rb = rb, pp = par(no.readonly = TRUE)), envir = .PlotPhyloEnv)
  }
  if (par.reset) {
    if (length(list(...))) {
      par(op)
    }
  }
  invisible(list(coords = ret$segs[-1, ], colorbreaks = colorbreaks, 
                 colordens = colorobj$colsdensity))
}


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
  
  #Switch back to plot.bammdata if the quantile(na.rm = TRUE) gets fixed in source
  xx = plot.bammdata2(edata, lwd=1, colorbreaks = colorbreaks, labels = T, lab.col = 'white')  
  cols = colorRampPalette(c('purple', 'white','springgreen'))(10)
  tiplabels(pch = 15, adj = (0.5 + offset), cex = 3, col = cols[home.regions$region])
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

burnInFrac = 0.2

sim4065 = output.unzip('../species-energy-simulation/raw_sim_output', 4065)
all.pops4065 = sim4065$all.populations
#extant.pops4065 = subset(all.pops4065, extant == 1)
#extant.phy4065 = read.tree('sim4065/extant_phy4065.tre')
#edata4065 = getEventData(extant.phy4065, eventdata = "sim4065/sim4065_event_data.txt", burnin = burnInFrac)
extant.pops4065.30k = subset(all.pops4065, time.of.origin <= 30000 & time.of.extinction > 30000)
extant.phy4065.30k = read.tree('sim4065-30k/run3/extant_phy4065_30k.tre')
edata4065.30k = getEventData(extant.phy4065.30k, eventdata = "sim4065-30k/run3/sim4065_30k.run3_event_data.txt", burnin = burnInFrac)

sim3465 = output.unzip('../species-energy-simulation/raw_sim_output', 3465)
all.pops3465 = sim3465$all.populations
extant.pops3465 = subset(all.pops3465, extant == 1)
extant.phy3465 = read.tree('sim3465/run1/extant_phy3465.tre')
edata3465 = getEventData(extant.phy3465, eventdata = "sim3465/run1/sim3465_event_data.txt", burnin = 0.05)

#sim3865 = output.unzip('../species-energy-simulation/raw_sim_output', 3865)
#all.pops3865 = sim3865$all.populations
#extant.pops3865 = subset(all.pops3865, extant == 1)
#extant.phy3865 = read.tree('sim3865/extant_phy3865.tre')
#edata3865 = getEventData(extant.phy3865, eventdata = "sim3865/sim3865_event_data.txt", burnin = burnInFrac)

sim5525 = output.unzip('../species-energy-simulation/raw_sim_output', 5525)
all.pops5525 = sim5525$all.populations
#extant.pops5525 = subset(all.pops5525, extant == 1)
#extant.phy5525 = read.tree('sim5525/extant_phy5525.tre')
#edata5525 = getEventData(extant.phy5525, eventdata = "sim5525/sim5525_event_data.txt", burnin = burnInFrac)
extant.pops5525.30k = subset(all.pops5525, time.of.origin <= 30000 & time.of.extinction > 30000)
extant.phy5525.30k = read.tree('sim5525-30k/run7/extant_phy5525_30k.tre')
edata5525.30k = getEventData(extant.phy5525.30k, eventdata = "sim5525-30k/run7/sim5525_30k.run7_event_data.txt", burnin = burnInFrac)

sim5625 = output.unzip('../species-energy-simulation/raw_sim_output', 5625)
all.pops5625 = sim5625$all.populations
#extant.pops5625 = subset(all.pops5625, extant == 1)
#extant.phy5625 = read.tree('sim5625/extant_phy5625.tre')
#edata5625 = getEventData(extant.phy5625, eventdata = "sim5625/sim5625_event_data.txt", burnin = burnInFrac)
extant.pops5625.30k = subset(all.pops5625, time.of.origin <= 30000 & time.of.extinction > 30000)
extant.phy5625.30k = read.tree('sim5625-30k/run2/extant_phy5625_30k.tre')
edata5625.30k = getEventData(extant.phy5625.30k, eventdata = "sim5625-30k/run2/sim5625_30k.run2_event_data.txt", burnin = burnInFrac)

# Use color legend of sim4065-30k for all plots
xx = plot.bammdata(edata4065.30k, show=F)

# Sims at t=30000, paired with rate through time plots
pdf(paste('rate_phylo_plot_4scenarios_', Sys.Date(), '.pdf', sep = ''),
    height = 10, width = 14)
par(mar = c(0, 4, 0, 2), oma = c(3, 3, 0, 0))
layout(matrix(c(rep(1:4, times = 4), 5:8), ncol=4, byrow=T))
ofs = 0.1

bamm.plot(extant.phy3465, edata3465, extant.pops3465, colorbreaks = xx$colorbreaks, offset = ofs)
mtext("a)", 2, las = 1, at = 1, outer=F, line = 1, cex = 1.5)

bamm.plot(extant.phy4065.30k, edata4065.30k, extant.pops4065.30k, offset = ofs)
mtext("b)", 2, las = 1, at = 1, outer=F, line = 1, cex = 1.5)

bamm.plot(extant.phy5525.30k, edata5525.30k, extant.pops5525.30k, colorbreaks = xx$colorbreaks, offset = ofs)
mtext("c)", 2, las = 1, at = 1, outer=F, line = 1, cex = 1.5)

bamm.plot(extant.phy5625.30k, edata5625.30k, extant.pops5625.30k, colorbreaks = xx$colorbreaks, offset = ofs)
mtext("d)", 2, las = 1, at = 1, outer=F, line = 1, cex = 1.5)

plotRateThroughTime(edata3465, axis.labels = F, yticks = 4, avgCol = 'black', intervalCol = 'black')

plotRateThroughTime(edata4065.30k, axis.labels = F, yticks = 4, avgCol = 'black', intervalCol = 'black')

plotRateThroughTime(edata5525.30k, axis.labels = F, yticks = 4, avgCol = 'black', intervalCol = 'black')

plotRateThroughTime(edata5625.30k, axis.labels = F, yticks = 4, avgCol = 'black', intervalCol = 'black')

mtext("Time steps before present", 1, outer=T, cex = 1.5)
mtext("Rate", 2, outer = T, adj = 0.12, cex = 1.5)

dev.off()









if(0) {

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
bamm.plot(extant.phy5625, edata5625, extant.pops5625)
mtext("Sim 5625, Disturbance gradient, t = 100k", 3)
#Page 4
bamm.plot(extant.phy3465, edata3465, extant.pops3465)
mtext("Sim 3465, Niche conservatism", 3)

dev.off()


####################################################
# Plotting tip-specific speciation rates versus region

plot.regLambda = function(simID, edata, extant.pops, fitline = NULL) {
  specnRates = data.frame(spp.name = edata$tip.label, Lambda = edata$meanTipLambda)
  extantRates = merge(extant.pops, specnRates, by = 'spp.name', all.x = T)
  regRates = aggregate(extantRates$Lambda, by = list(extantRates$region), function(x) mean(x, na.rm=T))
  plot(extantRates$region, extantRates$Lambda, xlab = '', ylab = '', xaxt = 'n')
  if (fitline == "spline") {
    points(smooth.spline(extantRates$region[is.na(extantRates$Lambda)==F], extantRates$Lambda[is.na(extantRates$Lambda)==F], df = 4),type='l',col='red')
  } else if (fitline == "linear") {
    abline(lm(extantRates$Lambda ~ extantRates$region), col = 'red')
  }
  mtext("Temperate", 1, adj = 0, cex = .75)
  mtext("Tropical", 1, adj = 1, cex = .75)
}

# t = 100k
pdf('specnRate_vs_latitude_4scenarios_100k.pdf', height = 8, width = 10)
par(mfrow = c(2, 2), oma = c(4, 4, 0, 0), mar = c(4, 5, 3, 1), las = 0)
plot.regLambda(3465, edata3465, extant.pops3465, "spline")
mtext("Niche conservatism", 3, line = 1)
plot.regLambda(4065, edata4065, extant.pops4065, "spline")
mtext("Energy gradient", 3, line = 1)
plot.regLambda(5525, edata5525, extant.pops5525, "spline")
mtext("Speciation gradient", 3, line = 1)
plot.regLambda(5625, edata5625, extant.pops5625, "spline")
mtext("Disturbance gradient", 3, line = 1)
mtext("Region", 1, outer=T, cex = 1.5)
mtext("Speciation rate", 2, outer=T, cex = 1.5)
mtext("t = 100000", 1, outer = T, at = .9, line = 3)
dev.off()

edata4065.30k.4 = getEventData(extant.phy4065.30k, eventdata = "sim4065-30k/run4/sim4065_30k.run4_event_data.txt", burnin = burnInFrac)
edata4065.30k.5 = getEventData(extant.phy4065.30k, eventdata = "sim4065-30k/run5/sim4065_30k.run5_event_data.txt", burnin = burnInFrac)
edata4065.30k.6 = getEventData(extant.phy4065.30k, eventdata = "sim4065-30k/run6/sim4065_30k.run6_event_data.txt", burnin = burnInFrac)

# t = 30k
pdf('specnRate_vs_latitude_4scenarios_30k.pdf', height = 8, width = 10)
par(mfrow = c(2, 2), oma = c(4, 4, 0, 0), mar = c(4, 5, 3, 1), las = 0)
plot.regLambda(3465, edata3465, extant.pops3465, "spline")
mtext("Niche conservatism", 3, line = 1)
plot.regLambda(4065, edata4065.30k, extant.pops4065.30k, "spline")
mtext("Energy gradient", 3, line = 1)
plot.regLambda(5525, edata5525.30k, extant.pops5525.30k, "spline")
mtext("Speciation gradient", 3, line = 1)
plot.regLambda(5625, edata5625.30k, extant.pops5625.30k, "spline")
mtext("Disturbance gradient", 3, line = 1)
mtext("Region", 1, outer=T, cex = 1.5)
mtext("Speciation rate", 2, outer=T, cex = 1.5)
mtext("t = 30000", 1, outer = T, at = .9, line = 3)
dev.off()

} #end if(0), unused code block