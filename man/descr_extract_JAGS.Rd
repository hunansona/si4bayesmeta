\name{descr_extract_JAGS}
\alias{descr_extract_JAGS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Extraction of descriptive statistics from objects provided by JAGS
%%  ~~function to do ... ~~
}
\description{Helper function for extraction of descriptive statistics (mean and sd) from output objects computed by \code{rjags}. See Roos et al. (2020) for details.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
descr_extract_JAGS(res_l, res_0, res_u)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{res_l}{a rjags object subject to a lower perturbation
%%     ~~Describe \code{res_l} here~~
}
  \item{res_0}{a base rjags object without any perturbation
%%     ~~Describe \code{res_0} here~~
}
  \item{res_u}{a rjags object subject to an upper perturbation 
%%     ~~Describe \code{res_u} here~~
}
}
\details{Extracts descriptive statistics (mean and sd) from \code{rjags} objects.
%%  ~~ If necessary, more details than the description above ~~
}
\value{A table with descriptive statistics (mean and sd) for \code{mu}, \code{log(tau)}, random effects \code{thetai} and \code{theta_new.}
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{Roos, M., Hunanyan, S., Bakka, H., Rue, H. (2020). Sensitivity and identification quantification by a relative latent model complexity perturbation in the Bayesian meta-analysis. Biometrical Journal (under revision).
%% ~put references to the literature/web site here ~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
