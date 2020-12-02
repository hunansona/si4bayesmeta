\name{rlmc_scaling_low}
\alias{rlmc_scaling_low}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{RLMC adjusted scaling of the study-specific within-study standard deviation in the likelihood: reduced impact of observations
%%  ~~function to do ... ~~
}
\description{Computes RLMC-adjusted scaled study-specific within-study standard deviations of observations in the likelihood according to Roos et al. (2020). The mean of the scaled observation is assumed to be fixed at the mean of the original observation and only its standard deviation is changed.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
rlmc_scaling_low(sigma, rlmc, hh)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sigma}{a vector of fixed standard deviations provided in the likelihood
%%     ~~Describe \code{sigma} here~~
}
  \item{rlmc}{the original RLMC
%%     ~~Describe \code{rlmc} here~~
}
  \item{hh}{perturbation for RLMC
%%     ~~Describe \code{hh} here~~
}
}
\details{In this function the reduced impact of observations by increasing their standard deviation is implemented, which we call the 'low' case. Larger scaled sigmai values mean that the observations are made less informative.
%%  ~~ If necessary, more details than the description above ~~
}
\value{Scaled sigma (a vector with within-study standard deviations).
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{Roos, M., Hunanyan, S., Bakka, H., Rue, H., (2020). Sensitivity and identification quantification by a relative latent model complexity perturbation in the Bayesian meta-analysis. Manuscript submitted to Biometrical Journal (under revision). 
%% ~put references to the literature/web site here ~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%%\examples{
%%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
