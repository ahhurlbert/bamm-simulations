# Function for grabbing summary stats regarding BAMM analyses based on the
# BAMM output in the specified folder

bamm_summary = function(simdir, edata = NULL, burnin = 0.2) {
  require(BAMMtools)
  require(coda)
  phyfile = list.files(simdir)[grep('phy', list.files(simdir))]
  phylo = read.tree(paste(simdir, '/', phyfile, sep = ''))
  mcfiles = list.files(simdir)[grep('mcmc_out', list.files(simdir))]
  mcmcout = read.csv(paste(simdir, '/', mcfiles[grep('^sim', mcfiles)], sep = ''))
  eventfile = list.files(simdir)[grep('event', list.files(simdir))]

  if (is.null(edata)) {
    edata = getEventData(phylo, eventdata = paste(simdir, '/', eventfile, sep = ''), 
                         burnin = burnin)
  }
    rtt = plotRateThroughTime(edata, plot = F)
  
  generations = max(mcmcout$generation) - 1
  writeFreq = mcmcout$generation[2] - mcmcout$generation[1]
  burnstart <- floor(burnin * nrow(mcmcout)) #burn in 
  postburn <- mcmcout[burnstart:nrow(mcmcout), ]
  effSize_Nshift = effectiveSize(postburn$N_shifts)
  effSize_logLik = effectiveSize(postburn$logLik)
  post_probs <- table(postburn$N_shifts) / nrow(postburn)
  modal_Nshift = as.numeric(names(which(post_probs == max(post_probs))))
  p_Nshift = max(post_probs)
  meanTipLambda = mean(edata$meanTipLambda)
  Lambda_max_over_min = max(rtt$avg)/min(rtt$avg)
  finalLogLik = mean(mcmcout$logLik[(nrow(mcmcout)-10):nrow(mcmcout)])
  
  output = data.frame(cbind(simdir, generations, writeFreq, burnin, effSize_Nshift, effSize_logLik,
                            modal_Nshift, p_Nshift, meanTipLambda, Lambda_max_over_min, finalLogLik))
  rm(list = c('edata', 'mcmcout'))
  return(output)
}

simdirs = c('sim3465', 
            'sim4065-30k', 
            'sim5525-30k', 
            'sim4065', 
            'sim3865', 
            'sim5525', 
            'sim5625', 
            'sim5625-30k')


summary.out = c()
for (s in simdirs) {
  tmp = bamm_summary(s)
  tmp$simID = s
  summary.out = rbind(summary.out, tmp)
}
summary = data.frame(summary.out)

