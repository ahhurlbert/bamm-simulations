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
plot.bammdata2 = function (x, method = "phylogram", vtheta = 5, rbf = 0.001, show = TRUE, 
                           labels = FALSE, legend = FALSE, spex = "s", lwd = 1, cex = 1, 
                           pal = "RdYlBu", colorbreaks = NULL, par.reset = TRUE, ...) 
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
                                                    ][phy$edge[, 2] == k, 1], adj = c(0, NA))
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
                                                      1, 4], phy$tip.label, cex = cex, pos = 4, offset = 0.25)
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
  xx = plot.bammdata2(edata, lwd=1, colorbreaks = colorbreaks)  
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


sim4065 = output.unzip('../species-energy-simulation/raw_sim_output', 4065)
sim3865 = output.unzip('../species-energy-simulation/raw_sim_output', 3865)
sim5525 = output.unzip('../species-energy-simulation/raw_sim_output', 5525)
sim3465 = output.unzip('../species-energy-simulation/raw_sim_output', 3465)

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


