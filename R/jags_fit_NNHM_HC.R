jags_fit_NNHM_HC<-function(df=df,
                           mu_prior_mean=0,
                           mu_prior_prec=1/(4^2),
                           tau_prior_prec=1/(0.5^2),
                           params=c("mu", "theta", "tau", "log_tau", "theta_new"),
                           nchains=4,
                           nadapt=4000,
                           nburnin=20000,
                           niter=120000,
                           nthin=4){
  
  # Computation of NNHM with a HC heterogeneity prior in JAGS
  
  
  jags_data=list(y=df$y, prec_y=1/df$sigma^2, 
                 mu_prior_mean=mu_prior_mean, mu_prior_prec=mu_prior_prec, 
                 tau_prior_prec=tau_prior_prec)
  
  # model 
  NNHM_HC_modelString <- "
# likelihood
model{
for(i in 1:length(y)){
y[i] ~ dnorm(theta[i], prec_y[i]);
theta[i] ~ dnorm(mu, tau_prec);
}

# predictive distribution 
theta_new ~ dnorm(mu, tau_prec);

# priors
mu ~ dnorm(mu_prior_mean, mu_prior_prec);
tau ~ dt(0, tau_prior_prec, 1)T(0,);
tau_prec<-1/(tau^2);
log_tau <- log(tau);
}
"

writeLines(NNHM_HC_modelString, con="TempModel.txt") # write to a file

# model initiation
NNHM_HC <- jags.model(
  file = "TempModel.txt", 
  data = jags_data,
  n.chains = nchains,
  n.adapt = nadapt
)


# burn-in
update(NNHM_HC, n.iter = nburnin)

# sampling/monitoring
NNHM_HC_j <- coda.samples(
  model = NNHM_HC, 
  variable.names = params, 
  n.iter = niter,
  thin = nthin
)


NNHM_HC_j_matrix<-as.matrix(NNHM_HC_j)

return(NNHM_HC_j_matrix)
}