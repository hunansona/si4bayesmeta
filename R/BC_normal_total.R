BC_normal_total <-
function(m1, sd1, m0, sd0){
  ## computation of the total BC under normality assumption BC_normal(sd_post_mu,sd_pri_mu)*BC_normal(mean_post_mu,mean_pri_mu) according to Roos et al. (2020)
  return(BC_normal_sd1sd0(sd1=sd1, sd0=sd0)*BC_normal_m1m0(m1=m1, sd1=sd1, m0=m0, sd0=sd0))
}
