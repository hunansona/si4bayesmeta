h_choice_all <-
function(df){
  ## Function for the choice of the numerical h step value for RLMC perturbations for actual data and for all 6
  ## scenarios given two heterogeneity priors HN and HC and three RLMC targets (0.25,0.5,0.75) in Roos et al. (2020)
  ## What typical RLMC perturbation does the epsilon grid correspond to for actual data and for all 6 scenarios?
  ## input
  ## df: data frame in bayesmeta format
  ### The choice of grid_epsilon: grid_epsilon=0.00354 adjusts for an earlier approach to epsilon-local sensitivity in Roos et al. (2015)
  ### Motivation for the choice of the grid_epsilon value in Roos et al. (2015)
  ### according to the calibration
  ### Nmuh1N012h(0.01) # 0.003535523
  ### epsilon for the local grid in Roos et al. (2015)
  ### grid_epsilon<-0.00354
  ### check the calibration
  ### h2Nmuh1N01(grid_epsilon) # 0.01001266
  ## output
  ## h: mean perturbation
  ## val: all perturbations
  
  # Get RLMC-adjusted scale parameters for HN and HC heterogeneity priors for RLMC targets (0.25,0.50,0.75)
  
  tau_HN_rlmc025_s<-pri_par_adjust_HN(df=df, rlmc=0.25, tail_prob=0.5)$p_HN
  tau_HN_rlmc050_s<-pri_par_adjust_HN(df=df, rlmc=0.50, tail_prob=0.5)$p_HN
  tau_HN_rlmc075_s<-pri_par_adjust_HN(df=df, rlmc=0.75, tail_prob=0.5)$p_HN
  
  tau_HC_rlmc025_s<-pri_par_adjust_HC(df=df, rlmc=0.25, tail_prob=0.5)$p_HC
  tau_HC_rlmc050_s<-pri_par_adjust_HC(df=df, rlmc=0.50, tail_prob=0.5)$p_HC
  tau_HC_rlmc075_s<-pri_par_adjust_HC(df=df, rlmc=0.75, tail_prob=0.5)$p_HC
  
  # settings
  grid_epsilon<-0.00354
  
  # Scenario 1: HN and RLMC=0.25
  # Compute epsilon grids for HN according to Ott et al. (2019)
  grid_tau_HN_rlmc025_s<-HN_A0_2_Al_Au(AA0=tau_HN_rlmc025_s, eps=grid_epsilon)
  # compute median RLMC values of the base parameter and both parameters in the grid
  r10<-median_rlmc(df=df, r.tau.prior=function(MM)rhalfnormal(n=MM, scale=tau_HN_rlmc025_s), MM=1000000)$median_rlmc_ref
  r1l<-median_rlmc(df=df, r.tau.prior=function(MM)rhalfnormal(n=MM, scale=grid_tau_HN_rlmc025_s[1]), MM=1000000)$median_rlmc_ref
  r1u<-median_rlmc(df=df, r.tau.prior=function(MM)rhalfnormal(n=MM, scale=grid_tau_HN_rlmc025_s[2]), MM=1000000)$median_rlmc_ref

  # Scenario 2: HN and RLMC=0.50
  # Compute epsilon grids for HN according to Ott et al. (2019)
  grid_tau_HN_rlmc050_s<-HN_A0_2_Al_Au(AA0=tau_HN_rlmc050_s, eps=grid_epsilon)
  # compute median RLMC values of the base parameter and both parameters in the grid
  r20<-median_rlmc(df=df, r.tau.prior=function(MM)rhalfnormal(n=MM, scale=tau_HN_rlmc050_s), MM=1000000)$median_rlmc_ref
  r2l<-median_rlmc(df=df, r.tau.prior=function(MM)rhalfnormal(n=MM, scale=grid_tau_HN_rlmc050_s[1]), MM=1000000)$median_rlmc_ref
  r2u<-median_rlmc(df=df, r.tau.prior=function(MM)rhalfnormal(n=MM, scale=grid_tau_HN_rlmc050_s[2]), MM=1000000)$median_rlmc_ref
  
  # Scenario 3: HN and RLMC=0.75
  # Compute epsilon grids for HN according to Ott et al. (2019)
  grid_tau_HN_rlmc075_s<-HN_A0_2_Al_Au(AA0=tau_HN_rlmc075_s, eps=grid_epsilon)
  # compute median RLMC values of the base parameter and both parameters in the grid
  r30<-median_rlmc(df=df, r.tau.prior=function(MM)rhalfnormal(n=MM, scale=tau_HN_rlmc075_s), MM=1000000)$median_rlmc_ref
  r3l<-median_rlmc(df=df, r.tau.prior=function(MM)rhalfnormal(n=MM, scale=grid_tau_HN_rlmc075_s[1]), MM=1000000)$median_rlmc_ref
  r3u<-median_rlmc(df=df, r.tau.prior=function(MM)rhalfnormal(n=MM, scale=grid_tau_HN_rlmc075_s[2]), MM=1000000)$median_rlmc_ref

  
  # Scenario 4: HC and RLMC=0.25
  # Compute epsilon grids for HN according to Ott et al. (2019)
  grid_tau_HC_rlmc025_s<-HC_A0_2_Al_Au(AA0=tau_HC_rlmc025_s, eps=grid_epsilon)
  # compute median RLMC values of the base parameter and both parameters in the grid
  r40<-median_rlmc(df=df, r.tau.prior=function(MM)bayesmeta::rhalfcauchy(n=MM, scale=tau_HC_rlmc025_s), MM=1000000)$median_rlmc_ref
  r4l<-median_rlmc(df=df, r.tau.prior=function(MM)bayesmeta::rhalfcauchy(n=MM, scale=grid_tau_HC_rlmc025_s[1]), MM=1000000)$median_rlmc_ref
  r4u<-median_rlmc(df=df, r.tau.prior=function(MM)bayesmeta::rhalfcauchy(n=MM, scale=grid_tau_HC_rlmc025_s[2]), MM=1000000)$median_rlmc_ref
  
  # Scenario 5: HC and RLMC=0.50
  # Compute epsilon grids for HN according to Ott et al. (2019)
  grid_tau_HC_rlmc050_s<-HC_A0_2_Al_Au(AA0=tau_HC_rlmc050_s, eps=grid_epsilon)
  # compute median RLMC values of the base parameter and both parameters in the grid
  r50<-median_rlmc(df=df, r.tau.prior=function(MM)bayesmeta::rhalfcauchy(n=MM, scale=tau_HC_rlmc050_s), MM=1000000)$median_rlmc_ref
  r5l<-median_rlmc(df=df, r.tau.prior=function(MM)bayesmeta::rhalfcauchy(n=MM, scale=grid_tau_HC_rlmc050_s[1]), MM=1000000)$median_rlmc_ref
  r5u<-median_rlmc(df=df, r.tau.prior=function(MM)bayesmeta::rhalfcauchy(n=MM, scale=grid_tau_HC_rlmc050_s[2]), MM=1000000)$median_rlmc_ref
  
  # Scenario 6: HC and RLMC=0.75
  # Compute epsilon grids for HN according to Ott et al. (2019)
  grid_tau_HC_rlmc075_s<-HC_A0_2_Al_Au(AA0=tau_HC_rlmc075_s, eps=grid_epsilon)
  # compute median RLMC values of the base parameter and both parameters in the grid
  r60<-median_rlmc(df=df, r.tau.prior=function(MM)bayesmeta::rhalfcauchy(n=MM, scale=tau_HC_rlmc075_s), MM=1000000)$median_rlmc_ref
  r6l<-median_rlmc(df=df, r.tau.prior=function(MM)bayesmeta::rhalfcauchy(n=MM, scale=grid_tau_HC_rlmc075_s[1]), MM=1000000)$median_rlmc_ref
  r6u<-median_rlmc(df=df, r.tau.prior=function(MM)bayesmeta::rhalfcauchy(n=MM, scale=grid_tau_HC_rlmc075_s[2]), MM=1000000)$median_rlmc_ref
  
  
  
  val<-c(abs(r1l-r10), abs(r1u-r10), abs(r2l-r20), abs(r2u-r20), abs(r3l-r30), abs(r3u-r30),
         abs(r4l-r40), abs(r4u-r40), abs(r5l-r50), abs(r5u-r50), abs(r6l-r60), abs(r6u-r60))
  h<-mean(val)

  return(list(h=h, val=val))
}
