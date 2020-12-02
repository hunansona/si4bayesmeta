sigma_ref <-
function(df){
  ## Function for computation of sigma_ref by geometric mean of sigmai
  ## computation of the reference standard deviation as suggested in equation (7) by Sorbye and Rue (2014)
  ## input:
  ## df: data frame with one column "sigma" containing the standard deviations sigmai in each study
  ## output:
  ## refernce standard deviation as suggested in equation (7) by Sorbye and Rue (2014)
  return(exp(mean(log(df$sigma))))
}
