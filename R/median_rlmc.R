median_rlmc <-
function(df,
                    r.tau.prior,
                    MM=1000000,
                    seed.value=12567){
  ## computation of the median RLMC by the MC simulation according to Ott et al. (2019)
  ## df: data frame containing a column df$sigma
  ## r.tau.prior: randomisation function for the prior
  ## MM: number of MC samples
  ## output:
  ## median_rlmc: MRLMC computed given individual study-specific sigmai values 
  ## median_rlmc_ref: MC MRLMC computed given a geometric mean (sigma_ref) of study-specific sigmai values
  ## Large discrepancy between median_rlmc and median_rlmc_ref can indicate outliers in study-specific sigmai values
  
  # supporting function
  sigma_ref<-function(df){
    ## sigma_ref by geometric mean of sigmai
    ## computation of the reference standard deviation as suggested in equation (7) by Sorbye and Rue (2014)
    ## input:
    ## df: data frame with one column "sigma" containing the standard deviations sigmai in each study
    ## output:
    ## refernce standard deviation as suggested in equation (7) by Sorbye and Rue (2014)
    return(exp(mean(log(df$sigma))))
  }
  
  # computations
  
  set.seed(seed.value)
  kk<-length(df$sigma) # check the number od studies in the data frame
  tau_sim<-r.tau.prior(MM) # generate a MC sample for tau
  
  # computations based on individual study-specific sigmai values
  pdsum<-0
  for (i in 1:kk){
    sim_ICCi<-tau_sim^2/(tau_sim^2+df$sigma[i]^2)
    pdsum<-pdsum+sim_ICCi
  }
  
  # computations based on sigma_ref values
  median_rlmc_ref<-median(tau_sim^2/(tau_sim^2+sigma_ref(df)^2))
  
  return(list(median_rlmc=median(pdsum/kk), median_rlmc_ref=median_rlmc_ref))
}
