\name{BC_normal_m1m0}
\alias{BC_normal_m1m0}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Bhattacharyya coefficient between two normal densities: the location modification 
%%  ~~function to do ... ~~
}
\description{Computes the Bhattacharyya coefficient (the mean-part) between two normal densities with means \code{m1} and \code{m0} and standard deviations \code{sd1} and \code{sd0} 
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
BC_normal_m1m0(m1, sd1, m0, sd0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{m1}{mean of the first normal density
%%     ~~Describe \code{m1} here~~
}
  \item{sd1}{standard deviation of the first normal density
%%     ~~Describe \code{sd1} here~~
}
  \item{m0}{mean of the second normal density
%%     ~~Describe \code{m0} here~~
}
  \item{sd0}{standard deviation of the second normal density
%%     ~~Describe \code{sd0} here~~
}
}
\details{The function returns the mean-part of the BC under normal approximation according to Roos et al. (2020). The mean-part is adjusted for standard deviations which corresponds to the Mahalanobis distance. This part quantifies the location modification.
%%  ~~ If necessary, more details than the description above ~~
}
\value{A numeric value is returned.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{Roos, M., Hunanyan, S., Bakka, H., Rue, H. (2020). Sensitivity and identification quantification by a relative latent model complexity perturbation in the Bayesian meta-analysis. Manuscript submitted to Biometrical Journal (under revision). 
%% ~put references to the literature/web site here ~
}

\examples{
BC_normal_m1m0(m1=0.4, sd1=0.3, m0=0.42, sd0=0.32)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
