h2Nmuh1N01 <-
function(hh){
  ## Normal calibration for a Hellinger distance value hh expressed in terms of muh for N(0,1) and N(muh,1) according to  Roos et al. (2015)
  return(sqrt(-8*log(1-hh^2)))
}
