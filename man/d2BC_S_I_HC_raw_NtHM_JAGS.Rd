\name{d2BC_S_I_HC_raw_NtHM_JAGS}
\alias{d2BC_S_I_HC_raw_NtHM_JAGS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Computation of Sensitivity and Identification measure expressed by the second derivative of the Bhattacharyya coefficient for a Bayesian NtHM with a HC heterogeneity prior in JAGS
%%  ~~function to do ... ~~
}
\description{Returns the second derivative of Bhattacharyya coefficient for a Bayesian NtHM with respect to Priori and Likelihood perturbations given a HC heterogeneity prior fitted by JAGS.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
d2BC_S_I_HC_raw_NtHM_JAGS(df, hh, rlmc = 0.5, mu_mean = 0, 
                              mu_sd = 4, nchains=4, nadapt=4000, 
                              nburnin=20000, niter=120000, 
                              nthin=4, tdf = 15)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{df}{original base data frame in a \code{bayesmeta} format
%%     ~~Describe \code{df} here~~
}
  \item{hh}{the step for the numerical computation of derivatives of BC with respect to changing RLMC in a Bayesian NtHM
%%     ~~Describe \code{hh} here~~
}
  \item{rlmc}{the value of the target RLMC, usually set to 0.25, 0.5 or 0.75
%%     ~~Describe \code{rlmc} here~~
}
  \item{mu_mean}{mean of the normal prior for mu
%%     ~~Describe \code{mu_mean} here~~
}
  \item{mu_sd}{sd of the normal prior for mu
%%     ~~Describe \code{mu_sd} here~~
}
\item{nchains}{number of parallel chains for the JAGS model. Default value is 4.
%%     ~~Describe \code{nchains} here~~
}
\item{nadapt}{number of iterations for adaptation. See \code{adapt} for details. If nadapt = 0 then no adaptation takes place. Default value is 4000.
%%     ~~Describe \code{nadapt} here~~
}
\item{nburnin}{number of burn-in iterations. Default value is 20000.
%%     ~~Describe \code{nburnin} here~~
}
\item{niter}{number of iterations to monitor. Default value is 120000.
%%     ~~Describe \code{niter} here~~
}
\item{nthin}{thinning interval for monitors. Default value is 4.
%%     ~~Describe \code{nburnin} here~~
}
\item{tdf}{degrees of freedom for Student t distribution. Default value is 15.
%%     ~~Describe \code{tdf} here~~
}
}
\details{Provides the Sensitivity and Identification estimates according to Roos et al. (2020). Sensitivity quantification is performed by a HC heterogeneity prior perturbation in a Bayesian NtHM with a fixed likelihood in JAGS. Identification quantification is obtained by likelihood perturbation in a Bayesian NtHM with a fixed HC heterogeneity prior. The measure is based on the computation of second derivatives of the Bhattacharyya coefficient with respect to RLMC evaluated at a base RLMC value. 
%%  ~~ If necessary, more details than the description above ~~
}
\value{Table with the S-I (sensitivity-identification) measure: with sensitivity ("S_d2BC_P") and identification ("I_d2BC_L").
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{Roos, M., Hunanyan, S., Bakka, H., Rue, H. (2020). Sensitivity and identification quantification by a relative latent model complexity perturbation in the Bayesian meta-analysis. Manuscript submitted to Biometrical Journal (under revision).
%% ~put references to the literature/web site here ~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{pri_par_adjust_HC}}, \code{\link{d2BC}}, \code{\link{descr_extract_JAGS}}, \code{\link{rjags}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
# Acute Graft rejection (AGR) data analyzed in Friede et al. (2017),  
# Sect. 3.2, URL: https://doi.org/10.1002/bimj.201500236
# First study: experimental group: 14 cases out of 61; 
# control group: 15 cases out of 20 
# Second study: experimental group: 4 cases out of 36; 
# control group: 11 cases out of 36  

rT<-c(14,4)
nT<-c(61,36)
rC<-c(15,11)
nC<-c(20,36)
  
df = data.frame(y = log((rT*(nC-rC))/(rC*(nT-rT))), 
            sigma = sqrt(1/rT+1/(nT-rT)+1/rC+1/(nC-rC)), 
                  labels = c(1:2))
  
# Set the normal prior parameter values for mu according to Roever
  #2018, bayesmeta
# Suitable for logOR and log-odds scale with the unit information 
  #argument  
mu_mean<-0
mu_sd<-4

#settings for JAGS modeling
nchains<-4
nadapt<-4000
nburnin<-20000
niter<-120000
nthin<-4
tdf <- 4
# set a global parameter for 
#RLMC perturbation

hh0<-0.0044 # value obtained empirically 
#by adjusting it for the 
#epsilon-local sensitivity approach
#with grid_epsilon=0.00354 for RLMC adjusted HN, HC 
#heterogeneity priors with targets (0.25,0.5,0.75) 
#for NtHM 

# Fix the value of the target RLMC
t_RLMC_025<-0.25

res_d2BC_S_I_rlmc025_HC_NtHM_JAGS_raw <- 
  d2BC_S_I_HC_raw_NtHM_JAGS(df=df, hh=hh0, rlmc=0.25, 
                              mu_mean=mu_mean, mu_sd=mu_sd,
                              nchains=nchains, nadapt=nadapt, 
                              nburnin=nburnin, niter=niter, 
                              nthin=nthin, tdf=tdf)

res_d2BC_S_I_rlmc025_HC_NtHM_JAGS_raw
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
