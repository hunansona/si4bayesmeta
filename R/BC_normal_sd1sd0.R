BC_normal_sd1sd0 <-
function(sd1, sd0){
  ## computation of BC_normal(sd_post_mu,sd_pri_mu) (the sd-part of the BC under normal approximation) according to Roos et al. (2020)
  ## this part quantifies the spread modification
  return(sqrt((2*sd1*sd0)/(sd1^2+sd0^2)))
}
