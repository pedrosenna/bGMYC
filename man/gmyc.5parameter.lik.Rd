\name{gmyc.5parameter.lik}
\alias{gmyc.5parameter.lik}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
5 parameter likelihood function}
\description{
returns the likelihood of the GMYC model if all 5 parameters are free to vary}
\usage{
gmyc.5parameter.lik(params, data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{params}{
a vector of five the five parameters in the model in the order y, c, py, pc, t}
  \item{data}{
output from gmyc.dataprep}
}
\details{
the Pons et al. implementation of the GMYC (and this one) treat the model as if it had 3 parameters, a yule rate change parameter, a coalescent rate change parameter and a threshold parameter. the yule and coalescent rates are fixed to an ML estimator at each possible threshold and not integrated over. allowing all 5 parameters to vary freely in some cases yields a flatter likelihood surface which in some cases (in simulation) has a peak very far away from the true value and very different from the 3 parameter implementation. using the ML estimator for the rate parameters yields empirically better peformance so we use it in our paper (Reid and Carstens 2012), but it took some time to figure that out. we include this function here for anyone who feels like messing with it.}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
