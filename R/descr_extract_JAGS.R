descr_extract_JAGS<-function (res_l, res_0, res_u){
  
  # extraction of descriptive statistics from objects provided by JAGS
  
  kk <- dim(res_0)[2]-4
  reff <- paste(rep("theta_", kk), c(1:kk), sep = "")
  names_row <- c("mu", "log_tau", reff, "theta_new")
  no_rows <- length(names_row)
  descr_names_col <- c("m_l", "sd_l", "m_0", 
                       "sd_0", "m_u", "sd_u")
  descr_collect <- matrix(NA, nrow = no_rows, ncol = length(descr_names_col), 
                          dimnames = list(names_row, descr_names_col))
  
  # res_l
  # descriptives mu
  descr_collect[1, 1] <- mean(res_l[, "mu"])
  descr_collect[1, 2] <- sd(res_l[, "mu"])
  
  # descriptives log_tau
  descr_collect[2, 1] <- mean(res_l[, "log_tau"])
  descr_collect[2, 2] <- sd(res_l[, "log_tau"])
  
  # descritives theta_i
  reff_name <- paste(rep("theta[", kk), c(1:kk), rep("]", kk), sep = "")
  for (i in 1:kk){
    descr_collect[i + 2, 1] <- mean(res_l[, reff_name[i]])
    descr_collect[i + 2, 2] <- sd(res_l[, reff_name[i]])
  }
  
  # descriptives theta_new
  descr_collect[kk + 3, 1] <- mean(res_l[, "theta_new"])
  descr_collect[kk + 3, 2] <- sd(res_l[, "theta_new"])
  
  
  # res_0
  # descriptives mu
  descr_collect[1, 3] <- mean(res_0[, "mu"])
  descr_collect[1, 4] <- sd(res_0[, "mu"])
  
  # descriptives log_tau
  descr_collect[2, 3] <- mean(res_0[, "log_tau"])
  descr_collect[2, 4] <- sd(res_0[, "log_tau"])
  
  # descritives theta_i
  reff_name <- paste(rep("theta[", kk), c(1:kk), rep("]", kk), sep = "")
  for (i in 1:kk){
    descr_collect[i + 2, 3] <- mean(res_0[, reff_name[i]])
    descr_collect[i + 2, 4] <- sd(res_0[, reff_name[i]])
  }
  
  # descriptives theta_new
  descr_collect[kk + 3, 3] <- mean(res_0[, "theta_new"])
  descr_collect[kk + 3, 4] <- sd(res_0[, "theta_new"])
  
  
  # res_u
  # descriptives mu
  descr_collect[1, 5] <- mean(res_u[, "mu"])
  descr_collect[1, 6] <- sd(res_u[, "mu"])
  
  # descriptives log_tau
  descr_collect[2, 5] <- mean(res_u[, "log_tau"])
  descr_collect[2, 6] <- sd(res_u[, "log_tau"])
  
  # descritives theta_i
  reff_name <- paste(rep("theta[", kk), c(1:kk), rep("]", kk), sep = "")
  for (i in 1:kk){
    descr_collect[i + 2, 5] <- mean(res_u[, reff_name[i]])
    descr_collect[i + 2, 6] <- sd(res_u[, reff_name[i]])
  }
  
  # descriptives theta_new
  descr_collect[kk + 3, 5] <- mean(res_u[, "theta_new"])
  descr_collect[kk + 3, 6] <- sd(res_u[, "theta_new"])
  
  
  return(descr_collect)
}
