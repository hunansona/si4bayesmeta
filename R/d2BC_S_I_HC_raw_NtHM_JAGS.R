d2BC_S_I_HC_raw_NtHM_JAGS<-function (df, hh, rlmc = 0.5, mu_mean = 0, mu_sd = 4,
                                     nchains=4, nadapt=4000, nburnin=20000,
                                     niter=120000, nthin=4, tdf=15){
  
  # Computation of S-I based on d2BC for NtHM with a HC heterogeneity prior in JAGS
  
  kk <- length(df$y)
  reff <- paste(rep("theta_", kk), c(1:kk), sep = "")
  names_row <- c("mu", "log_tau", reff, "theta_new")
  names_col_S <- c("S_d2BC_P")
  names_col_I <- c("I_d2BC_L")
  no_rows <- length(names_row)
  no_cols <- length(names_col_S)
  res_d2BC_P <- matrix(NA, nrow = no_rows, ncol = no_cols, 
                       dimnames = list(names_row, names_col_S))
  res_d2BC_L <- matrix(NA, nrow = no_rows, ncol = no_cols, 
                       dimnames = list(names_row, names_col_I))
  tau_rlmc_0 <- pri_par_adjust_HC(df = df, rlmc = rlmc, tail_prob = 0.5)$p_HC
  tau_rlmc_U <- pri_par_adjust_HC(df = df, rlmc = rlmc + hh, 
                                  tail_prob = 0.5)$p_HC
  tau_rlmc_L <- pri_par_adjust_HC(df = df, rlmc = rlmc - hh, 
                                  tail_prob = 0.5)$p_HC
  
  
  fit_P_0 <- jags_fit_NtHM_HC(df=df,
                              mu_prior_mean=mu_mean,
                              mu_prior_prec=1/(mu_sd^2),
                              tau_prior_prec=1/(tau_rlmc_0^2),
                              params=c("mu", "theta", "tau", "log_tau", "theta_new"),
                              nchains=nchains,
                              nadapt=nadapt,
                              nburnin=nburnin,
                              niter=niter,
                              nthin=nthin,
                              tdf=tdf)
  
  fit_L_0 <- fit_P_0
  
  fit_P_u <- jags_fit_NtHM_HC(df=df,
                              mu_prior_mean=mu_mean,
                              mu_prior_prec=1/(mu_sd^2),
                              tau_prior_prec=1/(tau_rlmc_U^2),
                              params=c("mu", "theta", "tau", "log_tau", "theta_new"),
                              nchains=nchains,
                              nadapt=nadapt,
                              nburnin=nburnin,
                              niter=niter,
                              nthin=nthin,
                              tdf=tdf)
  
  
  fit_P_l <- jags_fit_NtHM_HC(df=df,
                              mu_prior_mean=mu_mean,
                              mu_prior_prec=1/(mu_sd^2),
                              tau_prior_prec=1/(tau_rlmc_L^2),
                              params=c("mu", "theta", "tau", "log_tau", "theta_new"),
                              nchains=nchains,
                              nadapt=nadapt,
                              nburnin=nburnin,
                              niter=niter,
                              nthin=nthin,
                              tdf=tdf)
  
  
  res_d2BC_P[, 1] <- d2BC(descr_collect = descr_extract_JAGS(res_l = fit_P_l, 
                                                             res_0 = fit_P_0, res_u = fit_P_u), hh = hh)
  
  
  
  dfps_l <- data.frame(y = df$y, sigma = rlmc_scaling_low(sigma = df$sigma, 
                                                          rlmc = rlmc, hh = hh), labels = df$labels)
  dfps_u <- data.frame(y = df$y, sigma = rlmc_scaling_up(sigma = df$sigma, 
                                                         rlmc = rlmc, hh = hh), labels = df$labels)
  
  fit_L_l <- jags_fit_NtHM_HC(df=dfps_l,
                              mu_prior_mean=mu_mean,
                              mu_prior_prec=1/(mu_sd^2),
                              tau_prior_prec=1/(tau_rlmc_0^2),
                              params=c("mu", "theta", "tau", "log_tau", "theta_new"),
                              nchains=nchains,
                              nadapt=nadapt,
                              nburnin=nburnin,
                              niter=niter,
                              nthin=nthin,
                              tdf=tdf)
  
  
  fit_L_u <- jags_fit_NtHM_HC(df=dfps_u,
                              mu_prior_mean=mu_mean,
                              mu_prior_prec=1/(mu_sd^2),
                              tau_prior_prec=1/(tau_rlmc_0^2),
                              params=c("mu", "theta", "tau", "log_tau", "theta_new"),
                              nchains=nchains,
                              nadapt=nadapt,
                              nburnin=nburnin,
                              niter=niter,
                              nthin=nthin,
                              tdf=tdf)
  
  
  res_d2BC_L[, 1] <- d2BC(descr_collect = descr_extract_JAGS(res_l = fit_L_l, 
                                                             res_0 = fit_L_0, res_u = fit_L_u), hh = hh)
  return(-cbind(res_d2BC_P, res_d2BC_L))
}