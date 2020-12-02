\name{HN_A0_2_Al_Au}
\alias{HN_A0_2_Al_Au}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Epsilon-grid for a HN heterogeneity prior
%%  ~~function to do ... ~~
}
\description{Returns the epsilon-grid for HN(\code{A0}) heterogeneity prior according to Ott et al. (2020).
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
HN_A0_2_Al_Au(AA0, eps = grid_epsilon)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{AA0}{scale parameter of the base HN distribution
%%     ~~Describe \code{AA0} here~~
}
  \item{eps}{epsilon for the epsilon grid
%%     ~~Describe \code{eps} here~~
}
}
\details{Local epsilon-grid for the HN distribution which consists of two scale parameters \code{AAl} and \code{AAu} such that 
%%H(HC(AAl), HC(AA0))=eps and H(HC(AAu), HC(AA0))=eps,
\deqn{H(HN(AAl), HN(AA0))=eps,}
\deqn{H(HN(AAu), HN(AA0))=eps.}
See more details in Ott et al. (2020).
%%  ~~ If necessary, more details than the description above ~~
}
\value{A vector of two scale parameters \code{AAl} and \code{AAu}, that is, the local epsilon-grid for the HN distribution which consists of two scale parameters \code{AAl} and \code{AAu}.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{Ott, M., Hunanyan, S., Held, L., Roos, M. 2020. The relative latent model complexity adjustment for heterogeneity prior specification in Bayesian meta-analysis. Statistical Methods in Medical Research (in preparation).
%% ~put references to the literature/web site here ~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ \code{\link{pri_par_adjust_HC}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
%%\examples{
%%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
