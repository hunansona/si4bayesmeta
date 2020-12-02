\name{jags_fit_NNHM_HC}
\alias{jags_fit_NNHM_HC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Bayesian NNHM with a HC heterogeneity prior in JAGS
%%  ~~function to do ... ~~
}
\description{Fits a Bayesian NNHM with a RLMC-adjusted HC heterogeneity prior given a target RLMC according to Roos et al. (2020) in \code{rjags}.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
jags_fit_NNHM_HC(df=df, mu_prior_mean=0, mu_prior_prec=1/(4^2), 
                        tau_prior_prec=1/(0.5^2),
                        params=c("mu", "theta", "tau", "log_tau", "theta_new"),
                        nchains=4, nadapt=4000, nburnin=20000, niter=120000,
                        nthin=4)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{df}{original base data frame in a \code{rjags} format
%%     ~~Describe \code{df} here~~
}
  \item{mu_prior_mean}{mean of the normal prior for \code{mu}. Default value is 0.
%%     ~~Describe \code{df} here~~
}
  \item{mu_prior_prec}{precision of the normal prior for \code{mu}. Default value is 0.0625.
%%     ~~Describe \code{df} here~~
}
  \item{tau_prior_prec}{precision of the HC heterogeneity prior \code{tau}. Default value is 4.
%%     ~~Describe \code{df} here~~
}
  \item{params}{list of parameters to be estimated in NNHM. 
%%     ~~Describe \code{df} here~~
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
}
\value{A matrix with the results of a Bayesian NNHM with a RLMC-adjusted HC heterogeneity prior fitted in \code{rjags}.
%%  ~Describe the value returned 
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{Roos, M., Hunanyan, S., Bakka, H., Rue, H., (2020). Sensitivity and identification quantification by a relative latent model complexity perturbation in the Bayesian meta-analysis. Biometrical Journal (under revision). 
%% ~put references to the literature/web site here ~
}

\note{Needs to have the package \code{rjags} loaded.
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

