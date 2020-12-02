\name{Nmuh1N012h}
\alias{Nmuh1N012h}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Hellinger distance between normal and standard normal densities 
%%  ~~function to do ... ~~
}
\description{Hellinger distance between standard normal \eqn{N(0,1)} and unit variance normal \eqn{N(mu,1)} densities. \deqn{h(mu) = \sqrt{1- \exp{(-mu^2/8)}}.}
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Nmuh1N012h(mu)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mu}{mean shift
%%     ~~Describe \code{mu} here~~
}
}
\details{See more details in Roos et al., 2015, Appendix B4.
%%  ~~ If necessary, more details than the description above ~~
}
\value{Vector of Hellinger distances, consisting of real numbers in [0,1].
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

%%\examples{
%%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
