\name{median_rlmc}
\alias{median_rlmc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Median RLMC
%%  ~~function to do ... ~~
}
\description{Computes the median RLMC by the Monte Carlo simulation according to Ott et al. (2020)
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
median_rlmc(df, r.tau.prior, MM = 1e+06, seed.value = 12567)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{df}{data frame containing a column \code{sigma}
%%     ~~Describe \code{df} here~~
}
  \item{r.tau.prior}{randomisation function for the prior
%%     ~~Describe \code{r.tau.prior} here~~
}
  \item{MM}{number of Monte Carlo samples
%%     ~~Describe \code{MM} here~~
}
  \item{seed.value}{random seed for random number generator
%%     ~~Describe \code{seed.value} here~~
}
}
\details{Given the MC sample \eqn{{RLMC^{(j)}, j = 1, \ldots, J}}, the median RLMC is \deqn{MRLMC = Median({RLMC^{(j)}, j = 1, \ldots, J)}. } 
See more details in Ott et al. (2020).
%%  ~~ If necessary, more details than the description above ~~
}
\value{\item{median_rlmc}{MRLMC computed given individual study-specific sigmai values }
\item{median_rlmc_ref}{MC MRLMC computed given a geometric mean (sigma_ref) of study-specific sigmai values}
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{Ott, M., Hunanyan, S., Held, L., Roos, M. (2020). The relative latent model complexity adjustment for heterogeneity prior specification in Bayesian meta-analysis. Statistical Methods in Medical Research (in preparation). 
%% ~put references to the literature/web site here ~
}

\note{Large discrepancy between \code{median_rlmc} and \code{median_rlmc_ref} can indicate outliers in study-specific sigmai values.
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{sigma_ref}}
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
                  
median_rlmc(df = df, r.tau.prior = function(n) {rhalfnormal(n, scale = 0.5)})

median_rlmc(df = df, r.tau.prior = function(n) {rhalfcauchy(n, scale = 1)})

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
