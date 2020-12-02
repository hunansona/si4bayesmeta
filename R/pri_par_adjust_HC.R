pri_par_adjust_HC <-
function(df, rlmc=0.5, tail_prob=0.5){
  ## Function for an RLMC-based adjustment of the scale parameter for a HC distribution according to Ott et al. (2019)
  ## input:
  ## df: data frame
  ## rlmc: requested target RLMC
  ## output:
  ## RLMC-adjusted scale parameter for HN
  
  # supporting functions
  
  sigma_ref<-function(df){
    ## Function for computation of sigma_ref by geometric mean of sigmai
    ## computation of the reference standard deviation as suggested in equation (7) by Sorbye and Rue (2014)
    ## input:
    ## df: data frame with one column "sigma" containing the standard deviations sigmai in each study
    ## output:
    ## refernce standard deviation as suggested in equation (7) by Sorbye and Rue (2014)
    return(exp(mean(log(df$sigma))))
  }
  
  AA_from_Ualpha_HC <- function(UU, alpha){
    # parameter of the HC distribution given a threshold UU and tail probability alpha according to Ott et al. (2019)
    return(UU/tan(pi*(1-alpha)/2))
  }
  
  # computations
  tau_ref<-sqrt(rlmc/(1-rlmc))*sigma_ref(df)

  # parameter for HC
  p_HC <- AA_from_Ualpha_HC(tau_ref, alpha=tail_prob)
  
  return(list(p_HC=p_HC))
  
}
