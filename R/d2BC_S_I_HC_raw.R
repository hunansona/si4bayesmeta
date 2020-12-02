d2BC_S_I_HC_raw <-
function(df, hh, rlmc=0.5, 
                      mu_mean=0, mu_sd=4){
  ## d2BC wrt Priori and Likelihood perturbations given a HC heterogeneity prior: provides the raw Sensitivity-Identification measure
  ## according to Roos et al. (2020)
  ## HC heterogeneity prior
  ## Sensitivity quantification by Prior perturbation in a NNHM with a fixed Likelihood
  ## Identification quantification by Likelihood perturbation in a NNHM with a fixed Prior
  ## Computation of second derivatives of BC wrt to RLMC at base RLMC with normal approximation
  ## input: 
  ## df: original base data frame in a bayesmeta format
  ## hh: the step for the numerical computation of derivatives of BC with respect to changing relative model complexity in a NNHM
  ## rlmc: the value of the target relative latent model complexity (RLMC) Usually set to 0.25 or 0.5 (or 0.75)
  ## mu_mean: mean of the normal prior for mu
  ## mu_sd: sd of the normal prior for mu
  ## output: 
  ## Table with the S_I measure: with sensitivity ("S_d2BC_P") and identification ("I_d2BC_L")
  
  
  
  ## initialisation of matrices to store the results
  
  kk<-length(df$y)
  reff<-paste(rep("theta_",kk),c(1:kk), sep="")
  names_row<-c("mu", "log_tau", reff, "theta_new")
  names_col_S<-c("S_d2BC_P")
  names_col_I<-c("I_d2BC_L")
  no_rows<-length(names_row)
  no_cols<-length(names_col_S)
  
  # collect results in matrices for each P and L perturbation separately
  res_d2BC_P<-matrix(NA, nrow=no_rows, ncol=no_cols, 
                     dimnames=list(names_row,names_col_S))
  
  res_d2BC_L<-matrix(NA, nrow=no_rows, ncol=no_cols, 
                     dimnames=list(names_row,names_col_I))
  
  
  
  
  ## d2BC_P (Sensitivity) computation
  
  # prior perturbation
  
  # tau.prior: (base prior parametr values) function containing the function with the density of the prior for tau
  tau_rlmc_0<-pri_par_adjust_HC(df=df, rlmc=rlmc, tail_prob=0.5)$p_HC
  # 0.531027 (base HN parameter value for RLMC=0.25 for AGR)
  
  # tau.prior_U: more RLMC, imposes less smoothness and more overfitting, less shrinkage to 0, prior has more spread and put more mass at values away from 0
  ## more RLMC corresponds to a higher complexity ie. more spread in tau must be allowed
  tau_rlmc_U<-pri_par_adjust_HC(df=df, rlmc=rlmc+hh, tail_prob=0.5)$p_HC
  # 0.536125 (higher RLMC corresponds to a higher complexity ie. more spread in tau must be allowed for, there is less prior impact wrt. data, imposes less smoothness and more overfitting, less shrinkage to 0)
  
  # tau.prior_L: less RLMC, imposes more smoothness and less overfitting, more shrinkage to 0, prior has less spread and puts more mass close 0
  ## lower RLMC corresponds to a lower complexity ie. less spread in tau must be allowed
  tau_rlmc_L<-pri_par_adjust_HC(df=df, rlmc=rlmc-hh, tail_prob=0.5)$p_HC
  # 0.525929 (lower RLMC corresponds to a lower complexity ie. less spread in tau must be allowed for, more prior impact wrt. data, imposes more smoothness, more shrinkage to 0)
  
  
  # bayesmeta computation
  
  # base model
  # uses a tail-adjusted prior for tau
  fit_P_0 <- bayesmeta::bayesmeta(y=df[,"y"], 
                       sigma=df[,"sigma"],
                       labels=df[,"labels"],
                       mu.prior.mean=mu_mean, mu.prior.sd=mu_sd,
                       tau.prior=function(t){dhalfcauchy(t, scale=tau_rlmc_0)})
  
  # The base model without any L and P perturbation is the same
  fit_L_0<-fit_P_0
  
  # lower impact of the prior (SD of the HN prior is larger, wider) (in such a case the likelihood has more impact wrt prior)
  # less shrikage to 0 (weight on 0) is assumed, it corresponds to a larger relative latent model complexity
  # uses a tail-adjusted prior for tau
  fit_P_u <- bayesmeta::bayesmeta(y=df[,"y"], 
                       sigma=df[,"sigma"],
                       labels=df[,"labels"],
                       mu.prior.mean=mu_mean, mu.prior.sd=mu_sd,
                       tau.prior=function(t){dhalfcauchy(t, scale=tau_rlmc_U)})
  
  
  # higher impact of the prior (SD of the HN prior is smaller, narrower) (in such a case the likelihood has less impact wrt prior)
  # more shrinkage to 0 (more on 0) is assumed, it corresponds to a smaller relative latent model complexity
  # uses a tail-adjusted prior for tau
  fit_P_l <- bayesmeta::bayesmeta(y=df[,"y"], 
                       sigma=df[,"sigma"],
                       labels=df[,"labels"],
                       mu.prior.mean=mu_mean, mu.prior.sd=mu_sd,
                       tau.prior=function(t){dhalfcauchy(t, scale=tau_rlmc_L)})
  
  
  ### Sensitvity computation
  
  # raw Sensitivity measure expressed by d2BC
  res_d2BC_P[,1]<-d2BC(descr_collect=descr_extract(res_l=fit_P_l, res_0=fit_P_0, res_u=fit_P_u), hh=hh)
  

  
  
  
  
  
  ## d2BC_L (Identification) computation
  
  # df_l: df scaled with a lower impact of the likelihood (larger sigmais)
  # (lower weight of the likelihood) incerased standard deviation of observations sigmai
  # RLMC is being decreased
  # (precision of observations is decreased)
  dfps_l<-data.frame(y=df$y,
                     sigma=rlmc_scaling_low(sigma=df$sigma, rlmc=rlmc, hh=hh),
                     labels=df$labels)
  
  # df_u: df scaled with an upper (higher) impact of the likelihood (smaller sigmais)
  # (uppper, higher weight of the likelihood) decreased standard deviation of observations sigmai
  # RLMC is being increased
  # (precision of observations is increased)
  dfps_u<-data.frame(y=df$y,
                     sigma=rlmc_scaling_up(sigma=df$sigma, rlmc=rlmc, hh=hh),
                     labels=df$labels)
  
  
  
  # bayesmeta computation
  
  # This model has been already computed as a base model fit_P_0
  # # base model
  # # uses a tail-adjusted prior for tau
  # fit_L_0 <- bayesmeta(y=df[,"y"], 
  #                      sigma=df[,"sigma"],
  #                      labels=df[,"labels"],
  #                      mu.prior.mean=mu_mean, mu.prior.sd=mu_sd,
  #                      tau.prior=function(t){dhalfcauchy(t, scale=tau_rlmc_0)})
  
  # (lower impact of the likelihood, RLMC is being decreased)
  # fixed prior for tau
  fit_L_l <- bayesmeta::bayesmeta(y=dfps_l[,"y"], 
                       sigma=dfps_l[,"sigma"],
                       labels=dfps_l[,"labels"],
                       mu.prior.mean=mu_mean, mu.prior.sd=mu_sd,
                       tau.prior=function(t){dhalfcauchy(t, scale=tau_rlmc_0)})
  
  
  # (higher impact of the likelihood, RLMC is being increased)
  # fixed prior for tau
  fit_L_u <- bayesmeta::bayesmeta(y=dfps_u[,"y"], 
                       sigma=dfps_u[,"sigma"],
                       labels=dfps_u[,"labels"],
                       mu.prior.mean=mu_mean, mu.prior.sd=mu_sd,
                       tau.prior=function(t){dhalfcauchy(t, scale=tau_rlmc_0)})
  
  ### Identification computation
  
  # raw Identification measure expressed by d2BC
  res_d2BC_L[,1]<-d2BC(descr_collect=descr_extract(res_l=fit_L_l, res_0=fit_L_0, res_u=fit_L_u), hh=hh)
  

  
  return(-cbind(res_d2BC_P,res_d2BC_L))
  
}
