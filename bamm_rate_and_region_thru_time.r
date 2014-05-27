require(BAMMtools)

sim3865 = output.unzip('../species-energy-simulation/raw_sim_output', 3865)
timerich3865 = sim3865$time.richness

# Calculate the most temperate region occupied at any point in time
# (for calculating the tropical-most region for temperate origin sims,
# change 'min' to 'max')
minreg = sapply(unique(timerich3865$time), function(x) 
  min(timerich3865$region[timerich3865$time==x & 
                            timerich3865$spp.rich[timerich3865$time ==x] != 0])

extant.phy3865 = read.tree('sim3865/extant_phy3865.tre')
burnInFrac = 0.2
edata3865 = getEventData(extant.phy3865, eventdata = "sim3865/sim3865_event_data.txt", burnin = 0.05)


par(mar = c(5,5,1,5))
plotRateThroughTime(edata3865)
par(new=T)
plot(1:max(timerich3865$time), (10-minreg), type = 'n', xaxt = "n", yaxt = "n", xlab = "", ylab = "")
mtext(c("Tropical", "Temperate"), 4, adj = c(0.1, 0.9))
points(smooth.spline(1:30000, (10-minreg), df = 8), type = 'l')