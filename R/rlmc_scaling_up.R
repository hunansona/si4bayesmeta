rlmc_scaling_up <-
function(sigma, rlmc, hh){
  ## Function for RLMC adjusted scaling of the study-specific within-study standard deviation in the likelihood for NNHM according to Roos et al. (2020)
  ## up: incerased impact of observations (likelihood) by reducing their standard deviation
  ## function computes RLMC adjusted scaled standard deviations of observations in the likelihood
  ## Assumption: the mean of the scaled observation is fixed at the mean of the original one (only Standard deviations are changed)
  ## input:
  ## sigma: a vector of fixed standard deviations provided in the likelihood
  ## rlmc: the original RLMC
  ## hh: perturbation for RLMC
  ## output:
  ## su: scaled sigma (a vector with within-study standard deviations)
  RLMCu<-rlmc+hh # it increases model complexity
  fku<-sqrt((1-RLMCu)/RLMCu)/sqrt((1-rlmc)/rlmc) # the scaling factor should be smaller than 1 to provide more weight to data by making their sigma smaller
  # smaller sigma means that the observations are made more informative
  su<-sigma*fku
  return(su)
}
