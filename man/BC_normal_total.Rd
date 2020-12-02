\name{BC_normal_total}
\alias{BC_normal_total}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Bhattacharyya coefficient between two normal densities
%%  ~~function to do ... ~~
}
\description{Computes the total Bhattacharyya coefficient between two normal densities with means \code{m1} and \code{m0} and standard deviations \code{sd1} and \code{sd0} %%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
BC_normal_total(m1, sd1, m0, sd0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{m1}{mean of the fist normal density
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
\details{BC is a symmetric measure of affinity (Roos et al., 2015). 
\deqn{BC(\pi_b(\psi | y), \pi_a(\psi | y)) = \int_{-\infty}^{\infty} \sqrt{\pi_b(\psi | y) \pi_a(\psi | y)} d \psi.}
The total BC is the product of the mean part BC and the sd-part BC.
%%  ~~ If necessary, more details than the description above ~~
}
\value{A numeric value between 0 and 1. 
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{Roos, M., Martins, T., Held, L., Rue, H. (2015). Sensitivity analysis for Bayesian hierarchical models. Bayesian Analysis 10(2), 321-349. https://projecteuclid.org/euclid.ba/1422884977
%% ~put references to the literature/web site here ~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{BC_normal_m1m0}}, \code{\link{BC_normal_sd1sd0}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
BC_normal_m1m0(m1=0.4, sd1=0.3, m0=0.42, sd0=0.32)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
