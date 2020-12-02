descr_extract <-
function(res_l, res_0, res_u){
  ## Helper function for extraction of descriptive statistics from 3 output objects computed by bayesmeta in Roos et al. (2020)
  ## input:
  ## res_l, res_0, res_u: 3 output objects computed by bayesmeta
  ## output:
  ## descr_collect: a table with descriptive statistics for mu, log(tau), random effects thetai and theta_new
  
  kk<-length(res_0$theta[1,])
  reff<-paste(rep("theta_",kk),c(1:kk), sep="")
  names_row<-c("mu", "log_tau", reff, "theta_new")
  no_rows<-length(names_row)
  
  descr_names_col<-c("m_l","sd_l","m_0","sd_0","m_u","sd_u")
  descr_collect<-matrix(NA, nrow=no_rows, ncol=length(descr_names_col), 
                        dimnames=list(names_row,descr_names_col))
  
  # limit of integration
  # ilim<-Inf
  ilim<-710
  
  ### res_l: extraction of descriptive statistics
  descr_collect[1,c(1,2)]<-res_l$summary[c(3,4),2] # mu
  
  # descriptive statistics (under normality assumption) for perturbed posterior for log_tau
  mean_post_log_tau_l<-integrate(function(x){x*res_l$dposterior(tau=exp(x))*exp(x)}, lower = -ilim, upper = ilim)$value
  sd_post_log_tau_l<-sqrt(integrate(function(x){x^2*res_l$dposterior(tau=exp(x))*exp(x)}, lower = -ilim, upper = ilim)$value-mean_post_log_tau_l^2)
  descr_collect[2,c(1,2)]<-c(mean_post_log_tau_l,sd_post_log_tau_l) # log(tau)
  
  descr_collect[c(3:(kk+2)),c(1,2)]<-t(res_l$theta[c(5,6),]) # random effects theta_i
  descr_collect[kk+3,c(1,2)]<-res_l$summary[c(3,4),3] # theta_new
  
  
  ### res_0: extraction of descriptive statistics
  descr_collect[1,c(3,4)]<-res_0$summary[c(3,4),2] # mu
  
  # descriptive statistics (under normality assumption) for posterior for log_tau
  mean_post_log_tau_0<-integrate(function(x){x*res_0$dposterior(tau=exp(x))*exp(x)}, lower = -ilim, upper = ilim)$value
  sd_post_log_tau_0<-sqrt(integrate(function(x){x^2*res_0$dposterior(tau=exp(x))*exp(x)}, lower = -ilim, upper = ilim)$value-mean_post_log_tau_0^2)
  descr_collect[2,c(3,4)]<-c(mean_post_log_tau_0,sd_post_log_tau_0) # log(tau)
  
  descr_collect[c(3:(kk+2)),c(3,4)]<-t(res_0$theta[c(5,6),]) # random effects theta_i
  descr_collect[kk+3,c(3,4)]<-res_0$summary[c(3,4),3] # theta_new
  
  
  
  ### res_u: extraction of descriptive statistics
  descr_collect[1,c(5,6)]<-res_u$summary[c(3,4),2] # mu
  
  # descriptive statistics (under normality assumption) for perturbed posterior for log_tau
  mean_post_log_tau_u<-integrate(function(x){x*res_u$dposterior(tau=exp(x))*exp(x)}, lower = -ilim, upper = ilim)$value
  sd_post_log_tau_u<-sqrt(integrate(function(x){x^2*res_u$dposterior(tau=exp(x))*exp(x)}, lower = -ilim, upper = ilim)$value-mean_post_log_tau_u^2)
  descr_collect[2,c(5,6)]<-c(mean_post_log_tau_u,sd_post_log_tau_u) # log(tau)
  
  descr_collect[c(3:(kk+2)),c(5,6)]<-t(res_u$theta[c(5,6),]) # random effects theta_i
  descr_collect[kk+3,c(5,6)]<-res_u$summary[c(3,4),3] # theta_new
  
  
  return(descr_collect)
}
