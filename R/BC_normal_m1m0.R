BC_normal_m1m0 <-
function(m1, sd1, m0, sd0){
  ## computation of BC_normal(mean_post_mu,mean_pri_mu) (the mean-part of the BC under normal approximation) according to Roos et al. (2020)
  ## the mean-part is adjusted for standard deviations (it corresponds to the Mahalanobis distance)
  ## this part quantifies the location modification
  return(exp(-((m1-m0)^2/(4*(sd1^2+sd0^2)))))
}
