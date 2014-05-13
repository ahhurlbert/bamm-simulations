library(BAMMtools)
library(coda)
library(ape)

# BEFORE RUNNING BAMM #################################################

# Take simulation phylogeny and drop extinct taxa (this takes a long time for large trees)

sim = output.unzip('raw_sim_output', simID)
all.pops = sim$all.populations
extant.pops = subset(all.pops, extant == 1)
phy = sim$phylo.out

extant_phy = drop.tips(phy, phy$tip.label[!phy$tip.label %in% extant.pops$spp.name])
write.tree(extant_phy, 'pathname/extant_phy.tre')

# Get appropriate priors (a file is written to the working directory called 'myPriors.txt')
# These will then get pasted in the Priors block of the BAMM control file
setBAMMpriors(extant_phy)


# AFTER RUNNING BAMM ##################################################

mcmcout = read.csv('mcmc_out.txt', header = T)
plot(mcmcout$logLik ~ mcmcout$generation)

#evaluate above plot to determine fraction to discard in burn in
burnInFrac = 0.2
burnstart <- floor(burnInFrac * nrow(mcmcout)) #burn in 
postburn <- mcmcout[burnstart:nrow(mcmcout), ]

effectiveSize(postburn$N_shifts)
effectiveSize(postburn$logLik)

post_probs <- table(postburn$N_shifts) / nrow(postburn)
# Likelihood ratio for 1 shift vs 2 shifts:
post_probs[1] / post_probs[2]

tree <- read.tree("extant_phy3865.tre")
edata <- getEventData(tree, eventdata = "event_data.txt", burnin=0.2)
shift_probs <- summary(edata)

#Plot summary
plot.bammdata(edata, lwd=1)

# Here we will plot the 250th sample from the posterior:
index <- 250
e2 <- subsetEventData(edata, index = index)
plot.bammdata(e2, lwd=2)
addBAMMshifts(e2, cex=2)

#Clade node numbers (in ape format) at which shifts occurred
shiftnodes <- getShiftNodesFromIndex(edata, index = index)

#Plot rate through time for the entire phylogeny
plotRateThroughTime(edata, ratetype="speciation")

#Plot rate through time just for a subclade
plotRateThroughTime(edata, ratetype="speciation", node = shiftnodes)

# Credible shift set
shiftprior = read.table('shiftPrior_mcmc_out.txt', header=T, sep = ',')
shiftprior = shiftprior[burnstart:nrow(shiftprior), ]
priorshifts <- getBranchShiftPriors(tree, shiftprior)

css <- credibleShiftSet(edata, set.limit = 0.95, threshold = priorshifts)
css$number.distinct
summary(css)
plot.credibleshiftset(css)

# Marginal branch probabilities
marg_probs <- marginalShiftProbsTree(edata)
plot.phylo(marg_probs)