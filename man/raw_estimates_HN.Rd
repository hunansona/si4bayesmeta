\name{raw_estimates_HN}
\alias{raw_estimates_HN}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Raw estimates for the Bayesian meta-analysis expressed by a Bayesian NNHM with a HN heterogeneity prior in bayesmeta
%%  ~~function to do ... ~~
}
\description{Computes the raw estimates for a Bayesian NNHM with a RLMC-adjusted HN heterogeneity prior given a target RLMC according to Roos et al. (2020).
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
raw_estimates_HN(df, rlmc = 0.5, mu_mean = 0, mu_sd = 4)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{df}{original base data frame in a \code{bayesmeta} format
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

\note{Needs to have the package \code{bayesmeta} loaded.
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{pri_par_adjust_HN}}
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

# Fix the value of the target RMC
t_RLMC_025<-0.25

####---- HN: raw parameter estimates for the base NNHM ----####

raw_est_rlmc025_HN<-raw_estimates_HN(df=df, rlmc=t_RLMC_025, mu_mean=mu_mean, mu_sd=mu_sd)
raw_est_rlmc025_HN
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
