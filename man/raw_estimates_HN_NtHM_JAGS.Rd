\name{raw_estimates_HN_NtHM_JAGS}
\alias{raw_estimates_HN_NtHM_JAGS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Raw estimates of Bayesian NtHM with a HN heterogeneity prior in JAGS
%%  ~~function to do ... ~~
}
\description{Computes the raw posterior estimates for a Bayesian NtHM with a RLMC-adjusted HN heterogeneity prior given a target RLMC according to Roos et al. (2020).
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
raw_estimates_HN_NtHM_JAGS(df, rlmc = 0.5, mu_mean = 0, mu_sd = 4,
                                      nchains=4, nadapt=4000, nburnin=20000,
                                      niter=120000, nthin=4, tdf=15)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{df}{original base data frame in a \code{rjags} format
%%     ~~Describe \code{df} here~~
}
  \item{rlmc}{the value of the target RLMC, usually set to 0.25, 0.5 or 0.75. Default value is 0.5.
%%     ~~Describe \code{rlmc} here~~
}
  \item{mu_mean}{mean of the normal prior for \code{mu}. Default value is 0.
%%     ~~Describe \code{mu_mean} here~~
}
  \item{mu_sd}{standard deviation of the normal prior for \code{mu}. Default value is 4.
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
\value{A matrix with the raw estimates and shortest 95\% CrI and its length for \code{mu}, \code{tau}, \code{thetai}, \code{i=1,..,k} and \code{theta_new}.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{Roos, M., Hunanyan, S., Bakka, H., Rue, H., (2020). Sensitivity and identification quantification by a relative latent model complexity perturbation in the Bayesian meta-analysis. Manuscript submitted to Biometrical Journal (under revision). 
%% ~put references to the literature/web site here ~
}

\note{Needs to have the package \code{rjags} loaded.
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{pri_par_adjust_HN}}, \code{\link{jags_fit_NtHM_HN}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
rT<-c(14,4)
nT<-c(61,36)
rC<-c(15,11)
nC<-c(20,36)
  
df = data.frame(y = log((rT*(nC-rC))/(rC*(nT-rT))), # log-odds-ratio
                  sigma = sqrt(1/rT+1/(nT-rT)+1/rC+1/(nC-rC)), # SE(log-odds-ratio)
                  labels = c(1:2))
# Set the normal prior parameter values for mu according to Roever 2018 bayesmeta
# Suitable for logOR and log-odds scale with the unit information argument
mu_mean<-0
mu_sd<-4

#settings for JAGS modeling
nchains<-4
nadapt<-4000
nburnin<-20000
niter<-120000
nthin<-4
tdf<-4

# Fix the value of the target RLMC
t_RLMC_025<-0.25

####---- HN: raw parameter estimates for the base NtHM ----####
raw_est_rlmc025_HN_NtHM_JAGS <- raw_estimates_HN_NtHM_JAGS(df=df, rlmc=t_RLMC_025, 
                                mu_mean=mu_mean, mu_sd=mu_sd, nchains=nchains, 
                                nadapt=nadapt, nburnin=nburnin, niter=niter, 
                                nthin=nthin,tdf=tdf)
raw_est_rlmc025_HN_NtHM_JAGS
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
