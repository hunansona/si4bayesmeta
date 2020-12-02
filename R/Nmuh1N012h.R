Nmuh1N012h <-
function(mu){
  ## Hellinger distance between  N(0,1) and N(muh,1)
  return(sqrt(1-exp(-(mu^2)/8)))
}
