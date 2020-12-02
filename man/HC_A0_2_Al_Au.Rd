\name{HC_A0_2_Al_Au}
\alias{HC_A0_2_Al_Au}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Epsilon-grid for a HC heterogeneity prior
%%  ~~function to do ... ~~
}
\description{Returns the epsilon-grid for HC(\code{A0}) heterogeneity prior according to Ott et al. (2020).
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
HC_A0_2_Al_Au(AA0, eps = grid_epsilon)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{AA0}{scale parameter of the base HC distribution
%%     ~~Describe \code{AA0} here~~
}
  \item{eps}{epsilon for the epsilon grid
%%     ~~Describe \code{eps} here~~
}
}
\details{Local epsilon-grid for the HC distribution which consists of two scale parameters \code{AAl} and \code{AAu} such that 
%%H(HC(AAl), HC(AA0))=eps and H(HC(AAu), HC(AA0))=eps,
\deqn{H(HC(AAl), HC(AA0))=eps,}
\deqn{H(HC(AAu), HC(AA0))=eps.}
See more details in Ott et al. (2020).
%%  ~~ If necessary, more details than the description above ~~
}
\value{A vector of two scale parameters \code{AAl} and \code{AAu}, that is, the local epsilon-grid for the HC distribution which consists of two scale parameters \code{AAl} and \code{AAu}.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{Ott, M., Hunanyan, S., Held, L., Roos, M. (2020). The relative latent model complexity adjustment for heterogeneity prior specification in Bayesian meta-analysis. Statistical Methods in Medical Research (in preparation). 
%% ~put references to the literature/web site here ~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{pri_par_adjust_HC}}
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
  
tau_HC_rlmc025_s <- pri_par_adjust_HC(df=df, rlmc=0.25, tail_prob=0.5)$p_HC
HC_A0_2_Al_Au(AA0=tau_HC_rlmc025_s, eps=0.00354)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
