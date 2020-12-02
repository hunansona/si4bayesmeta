raw_estimates_HN <-
function(df, rlmc=0.5, mu_mean=0, mu_sd=4){
  ## raw estimates for Bayesian NNHM meta-analyses with a HN heterogeneity prior subject to a RLMC-based adjustment
  ## given a target RLMC in Roos et al. (2020)
  ## input: 
  ## df: original base data frame in a bayesmeta format
  ## rlmc: the value of the target relative latent model complexity (RLMC) Usually set to 0.25 or 0.5 (or 0.75)
  ## mu_mean: mean of the normal prior for mu
  ## mu_sd: sd of the normal prior for mu
  ## output: 
  ## raw_estimates: a table with the raw estimates and shortest 95%CrI and ist length for mu, tau, thetai, i=1,..,k and theta_new

  
  ## initialisation of matrices to store the results
  
  kk<-length(df$y)
  reff<-paste(rep("theta_",kk),c(1:kk), sep="")
  names_row<-c("mu", "tau", reff, "theta_new")
  names_col<-c("estimate_mean", "estimate_sd", "estimate_CI_low", "estimate_CI_up", "length_estimate_CI")
  no_rows<-length(names_row)
  no_cols<-length(names_col)
  
  # collect results in a matrix for a base/target RLMC
  raw_estimates<-matrix(NA, nrow=no_rows, ncol=no_cols, 
                     dimnames=list(names_row,names_col))
  
  
  # RLMC-based adjustment
  
  # tau.prior: (base prior parametr values) function containing the function with the density of the prior for tau
  tau_rlmc_0<-pri_par_adjust_HN(df=df, rlmc=rlmc, tail_prob=0.5)$p_HN
  
  # bayesmeta computation
  
  # base model
  # uses a RLMC-based adjustment for the heterogeneity prior for tau
  res_0 <- bayesmeta::bayesmeta(y=df[,"y"], 
                       sigma=df[,"sigma"],
                       labels=df[,"labels"],
                       mu.prior.mean=mu_mean, mu.prior.sd=mu_sd,
                       tau.prior=function(t){bayesmeta::dhalfnormal(t, scale=tau_rlmc_0)})
  
  ### res_0: extraction of raw posterior estimates
  raw_estimates[1,c(1:4)]<-res_0$summary[c(3:6),2] # mu
  raw_estimates[2,c(1:4)]<-res_0$summary[c(3:6),1] # tau
  raw_estimates[c(3:(kk+2)),c(1:4)]<-t(res_0$theta[c(5:8),]) # random effects theta_i
  raw_estimates[kk+3,c(1:4)]<-res_0$summary[c(3:6),3] # theta_new
  
  # computation of the length of the credible interval
  raw_estimates[,5]<-raw_estimates[,"estimate_CI_up"]-raw_estimates[,"estimate_CI_low"]
  
  return(raw_estimates)
}
