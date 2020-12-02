HC_A0_2_Al_Au <-
function(AA0, eps=grid_epsilon){
  grid_epsilon <- NULL
  ## epsilon-grid for HC(A0) according to Ott et al. (2019)
  ## AA0: scale parameter of the base half Cauchy distribution
  ## eps: epsilon for the epsilon local grid
  ## output AAl, AAu: local epsilon grid for the half Cauchy distribution which consists of two scale parameters AAl and AAu
  ## such that H(HC(AAl), HC(AA0))=eps and H(HC(AAu), HC(AA0))=eps
  
  # HC
  obj_HC <- function(x, AA0, eps=eps){
    # Function for numerical search for the epsilon local grid for the scaled half cauchy distributions
    #library(bayesmeta)
    return(integrate(function(t) {sqrt(dhalfcauchy(t, scale=x)*dhalfcauchy(t, scale=AA0))}, lower = 0, upper = Inf)$value-(1-eps^2))
  }
  
  AAl <- uniroot(obj_HC, lower=0.0001, upper=AA0, tol= 1e-9, AA0=AA0, eps=eps)$root
  AAu <- uniroot(obj_HC, lower=AA0, upper=100, tol= 1e-9, AA0=AA0, eps=eps)$root
  
  return(c(AAl,AAu))
}
