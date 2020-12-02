\name{h2Nmuh1N01}
\alias{h2Nmuh1N01}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Calibration of the Hellinger distance
%%  ~~function to do ... ~~
}
\description{Returns the calibration of the Hellinger distance \code{hh} expressed as a shift in the mean between a standard normal distribution and a unit variance normal distribution with this mean.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
h2Nmuh1N01(hh)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{hh}{Hellinger distance value (or a vector of Hellinger distance values), a real value between 0 and 1.
%%     ~~Describe \code{hh} here~~
}
}
\details{Returns the calibration 
\deqn{\mu(hh) = \sqrt{-8 \log {(1-hh^2)}}}
for the given Hellinger distance \code{hh}, such that 
\deqn{H(N(\mu(hh), 1), N(0,1)) = hh.}
See Roos et al. (2015), Sect. 2.2 for more details.
%%  ~~ If necessary, more details than the description above ~~
}
\value{A mean value (or a vector of mean values). 
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

\examples{
h2Nmuh1N01(h=c(0.1, 0.5, 0.9))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
