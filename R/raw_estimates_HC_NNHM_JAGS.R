raw_estimates_HC_NNHM_JAGS<-function (df, rlmc = 0.5, mu_mean = 0, mu_sd = 4,
                                      nchains=4, nadapt=4000, nburnin=20000,
                                      niter=120000, nthin=4){
  
  # Computation of posterior estimates of NNHM with a HC heterogeneity prior in JAGS
  
  
  kk <- length(df$y)
  reff <- paste(rep("theta_", kk), c(1:kk), sep = "")
  names_row <- c("mu", "tau", reff, "theta_new")
  names_col <- c("estimate_mean", "estimate_sd", 
                 "estimate_CI_low", "estimate_CI_up", "length_estimate_CI")
  no_rows <- length(names_row)
  no_cols <- length(names_col)
  raw_estimates <- matrix(NA, nrow = no_rows, ncol = no_cols, 
                          dimnames = list(names_row, names_col))
  tau_rlmc_0 <- pri_par_adjust_HC(df = df, rlmc = rlmc, tail_prob = 0.5)$p_HC
  
  res_0 <- jags_fit_NNHM_HC(df=df,
                            mu_prior_mean=mu_mean,
                            mu_prior_prec=1/(mu_sd^2),
                            tau_prior_prec=1/(tau_rlmc_0^2),
                            params=c("mu", "theta", "tau", "log_tau", "theta_new"),
                            nchains=nchains,
                            nadapt=nadapt,
                            nburnin=nburnin,
                            niter=niter,
                            nthin=nthin)
  
  # descriptives mu
  raw_estimates[1, 1] <- mean(res_0[, "mu"])
  raw_estimates[1, 2] <- sd(res_0[, "mu"])
  raw_estimates[1, 3] <- quantile(res_0[, "mu"], probs=0.025)
  raw_estimates[1, 4] <- quantile(res_0[, "mu"], probs=0.975)
  
  # descriptives tau
  raw_estimates[2, 1] <- mean(res_0[, "tau"])
  raw_estimates[2, 2] <- sd(res_0[, "tau"])
  raw_estimates[2, 3] <- quantile(res_0[, "tau"], probs=0.025)
  raw_estimates[2, 4] <- quantile(res_0[, "tau"], probs=0.975)
  
  # descritives theta_i
  reff_name <- paste(rep("theta[", kk), c(1:kk), rep("]", kk), sep = "")
  for (i in 1:kk){
    raw_estimates[i + 2, 1] <- mean(res_0[, reff_name[i]])
    raw_estimates[i + 2, 2] <- sd(res_0[, reff_name[i]])
    raw_estimates[i + 2, 3] <- quantile(res_0[, reff_name[i]], probs=0.025)
    raw_estimates[i + 2, 4] <- quantile(res_0[, reff_name[i]], probs=0.975)
    
  }
  
  # descriptives theta_new
  raw_estimates[kk + 3, 1] <- mean(res_0[, "theta_new"])
  raw_estimates[kk + 3, 2] <- sd(res_0[, "theta_new"])
  raw_estimates[kk + 3, 3] <- quantile(res_0[, "theta_new"], probs=0.025)
  raw_estimates[kk + 3, 4] <- quantile(res_0[, "theta_new"], probs=0.975)
  
  raw_estimates[, 5] <- raw_estimates[, "estimate_CI_up"] - 
    raw_estimates[, "estimate_CI_low"]
  return(raw_estimates)
}