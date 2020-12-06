# si4bayesmeta
Sensitivity and identification for the Bayesian Meta-Analysis

Sensitivity and identification estimates for all the parameters in the Bayesian normal-normal hierarchical model (NNHM) and in the Bayesian normal-Student t hierarchical model (NtHM) induced by a Half-Normal (HN) and a Half-Cauchy (HC) heterogeneity priors. Estimates are produced by 4 functions for Bayesian NNHM fitted either with bayesmeta or JAGS and by 2 functions for Bayesian NtHM fitted with JAGS. Six scenarios are considered: target relative latent model complexity (RLMC) values fixed at 0.25, 0.5 and 0.75 with RLMC-adjusted HN and HC heterogeneity priors. Corresponding posterior estimates can be obtained by the functions of raw_estimates_type. The methodology implemented in this package has been developed in Roos et al. (2020).


*********************************************************************************
How to install an R package that is sitting on github to R (for example, from Rstudio)?

> install.packages("devtools")
> library("devtools")
> install_github("hunansona/si4bayesmeta") 

You might need 
> Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS="true")
> install_github("hunansona/si4bayesmeta") 


