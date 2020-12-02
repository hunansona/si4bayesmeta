\name{pri_par_adjust_HC}
\alias{pri_par_adjust_HC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{RLMC-based adjustment of the scale parameter for a HC distribution
%%  ~~function to do ... ~~
}
\description{Performs a RLMC-based adjustment of the scale parameter for a HC distribution according to Ott et al. (2020). 
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
pri_par_adjust_HC(df, rlmc = 0.5, tail_prob = 0.5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{df}{data frame in \code{bayesmeta} format
%%     ~~Describe \code{df} here~~
}
  \item{rlmc}{requested target RLMC
%%     ~~Describe \code{rlmc} here~~
}
  \item{tail_prob}{tail probability of a heterogeneity prior
%%     ~~Describe \code{tail_prob} here~~
}
}
\details{Given the target amount of RLMC, we use a reference threshold \eqn{U_{ref}} for heterogeneity prior adjustment
\deqn{U_{ref} = \sigma_{ref} \sqrt{RLMC/(1-RLMC)},}
where \eqn{\sigma_{ref}} is the reference standard deviation. 
The RLMC-adjustment on the median of different heterogeneity priors aligned at \eqn{U_{ref}} reference threshold
\deqn{P[\tau > U_{ref}] = \alpha = 50 \% . }
See more details in Ott et al. (2020).
%%  ~~ If necessary, more details than the description above ~~
}
\value{A list of
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{p_HC}{RLMC-based adjusted scale parameter for a HC distribution}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{Ott, M., Hunanyan, S., Held, L., Roos, M. (2020). The relative latent model complexity adjustment for heterogeneity prior specification in Bayesian meta-analysis. Statistical Methods in Medical Research (in preparation).
%% ~put references to the literature/web site here ~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{sigma_ref}}, package \code{\link{bayesmeta}}
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

tau_HC_rlmc025_s<-pri_par_adjust_HC(df=df, rlmc=0.25, tail_prob=0.5)
tau_HC_rlmc025_s
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
