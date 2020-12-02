\name{h_choice_all}
\alias{h_choice_all}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{The choice of the numerical \code{h} step value for RLMC perturbations
%%  ~~function to do ... ~~
}
\description{Function for the choice of the numerical \code{h} step value for RLMC perturbations for actual data and for all 6 scenarios given two RLMC-adjusted heterogeneity priors HN and HC and three RLMC targets (0.25, 0.5, 0.75) in Roos et al. (2020).
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
h_choice_all(df)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{df}{data frame in bayesmeta format
%%     ~~Describe \code{df} here~~
}
}
\details{See the motivation for the choice of the grid_epsilon value (\code{grid_epsilon}=0.00354) in Roos et al. (2015)
%%  ~~ If necessary, more details than the description above ~~
}
\value{
 \item{h}{mean perturbation}
  \item{val}{all perturbations}
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{Roos, M., Hunanyan, S., Bakka, H., Rue, H. (2020). Sensitivity and identification quantification by a relative latent model complexity perturbation in the Bayesian meta-analysis. Manuscript submitted to Biometrical Journal (under revision). 

Roos, M., Martins, T., Held, L., Rue, H. (2015). Sensitivity analysis for Bayesian hierarchical models. Bayesian Analysis 10(2), 321-349. https://projecteuclid.org/euclid.ba/1422884977
%% ~put references to the literature/web site here ~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{h2Nmuh1N01}}
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

h_choice_all(df=df)$h # AGR: 0.0044249
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
