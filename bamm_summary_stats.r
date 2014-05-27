
#NOTE: this function doesn't currently work for sims like 4065.30k because of difference between 
# directory name and file names

bamm_summary = function(mcmcout, phylo, simID, burnin = 0.2) {
  require(BAMMtools)
  require(coda)
  edata = getEventData(phylo, eventdata = paste('sim', simID, '/sim', simID, '_event_data.txt', sep = ''), burnin = burnin)
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
  
  output = data.frame(cbind(generations, writeFreq, burnin, effSize_Nshift, effSize_logLik,
                            modal_Nshift, p_Nshift, meanTipLambda, Lambda_max_over_min))
  return(output)
  
}
