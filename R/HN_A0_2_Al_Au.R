HN_A0_2_Al_Au <-
function(AA0, eps=grid_epsilon){
  grid_epsilon <- NULL
  ## epsilon-grid for HN(A0) according to Ott et al. (2019)
  ## AA0: scale parameter of the base half normal distribution
  ## eps: epsilon for the epsilon local grid
  ## output AAl, AAu: local epsilon grid for the half normal distribution which consists of two scale parameters AAl and AAu
  ## such that H(HN(AAl), HN(AA0))=eps and H(HN(AAu), HN(AA0))=eps
  
  AAl <-AA0*(1/(1-eps^2)^2-sqrt(1/(1-eps^2)^4-1))
  AAu<-AA0*(1/(1-eps^2)^2+sqrt(1/(1-eps^2)^4-1))
  return(c(AAl,AAu))
}
