rlmc_scaling_low <-
function(sigma, rlmc, hh){
  ## Function for RLMC adjusted scaling of the study-specific within-study standard deviation in the likelihood for NNHM according to Roos et al. (2020)
  ## low: reduced impact of observations by increasing their standard deviation
  ## function computes RLMC adjusted scaled standard deviations of observations in the likelihood
  ## Assumption: the mean of the scaled observation is fixed at the mean of the original one (only Standard deviations are changed)
  ## input:
  ## sigma: a vector of fixed standard deviations provided in the likelihood
  ## rlmc: the original RLMC
  ## hh: perturbation for RLMC
  ## output:
  ## sl: scaled sigma (a vector with within-study standard deviations)
  RLMCl<-rlmc-hh # it decreases model complexity
  fkl<-sqrt((1-RLMCl)/RLMCl)/sqrt((1-rlmc)/rlmc) # scaling factor should be larger than 1 to provide less weight to data by making their sigma larger
  # larger sigma means that the observations are made less informative
  sl<-sigma*fkl
  return(sl)
}
