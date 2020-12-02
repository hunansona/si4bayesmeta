jags_fit_NtHM_HN<-function(df=df,
                           mu_prior_mean=0,
                           mu_prior_prec=1/(4^2),
                           tau_prior_prec=1/(0.5^2),
                           params=c("mu", "theta", "tau", "log_tau", "theta_new"),
                           nchains=4,
                           nadapt=4000,
                           nburnin=20000,
                           niter=120000,
                           nthin=4,
                           tdf=15){
  
  # Computation of NtHM with a HN heterogeneity prior in JAGS 
  
  jags_data=list(y=df$y, prec_y=1/df$sigma^2, 
                 mu_prior_mean=mu_prior_mean, mu_prior_prec=mu_prior_prec, 
                 tau_prior_prec=tau_prior_prec, tdf=tdf)
  
  # model 
  NtHM_HN_modelString <- "
# likelihood
model{
for(i in 1:length(y)){
y[i] ~ dnorm(theta[i], prec_y[i]);
theta[i] ~ dt(mu, tau_prec, tdf);
}

# predictive distribution 
theta_new ~ dnorm(mu, tau_prec);

# priors
mu ~ dnorm(mu_prior_mean, mu_prior_prec);
tau ~ dnorm(0, tau_prior_prec)T(0,);
tau_prec<-1/(tau^2);
log_tau <- log(tau);
}
"

writeLines(NtHM_HN_modelString, con="TempModel.txt") # write to a file

# model initiation
NtHM_HN <- jags.model(
  file = "TempModel.txt", 
  data = jags_data,
  n.chains = nchains,
  n.adapt = nadapt
)


# burn-in
update(NtHM_HN, n.iter = nburnin)

# sampling/monitoring
NtHM_HN_j <- coda.samples(
  model = NtHM_HN, 
  variable.names = params, 
  n.iter = niter,
  thin = nthin
)


NtHM_HN_j_matrix<-as.matrix(NtHM_HN_j)

return(NtHM_HN_j_matrix)
}
