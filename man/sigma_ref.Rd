\name{sigma_ref}
\alias{sigma_ref}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Reference within-study standard deviation
%%  ~~function to do ... ~~
}
\description{Computes the reference within-study standard deviation through the geometric mean as suggested in equation (7) by Sorbye and Rue (2014).
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
sigma_ref(df)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{df}{data frame with one column "sigma" containing the standard deviations sigmai in each study
%%     ~~Describe \code{df} here~~
}
}

\value{Refernce standard deviation expressed by the geometric mean as suggested in equation (7) by Sorbye and Rue (2014). Non-negative real number.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{Sorbye, S., Rue, H. (2014). Scaling intrinsic Gaussian Markov random field priors in spatial modelling. Spatial Statistics 8, 39-51. https://doi.org/10.1016/j.spasta.2013.06.004
%% ~put references to the literature/web site here ~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
rT<-c(14,4)
nT<-c(61,36)
rC<-c(15,11)
nC<-c(20,36)
  
df = data.frame(y = log((rT*(nC-rC))/(rC*(nT-rT))), # log-odds-ratio
                sigma = sqrt(1/rT+1/(nT-rT)+1/rC+1/(nC-rC)), # SE(log-odds-ratio)
                labels = c(1:2))
  
sigma_ref(df) # 0.620372
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
