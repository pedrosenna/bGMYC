\name{bgmyc.lik}
\alias{bgmyc.lik}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
gmyc likelihood function}
\description{
returns the log likelihood of the model given the data and parameter values}
\usage{
bgmyc.lik(params, data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{params}{
a vector of parameter values to be evaluated in the order: py, pc, t}
  \item{data}{
an output from gmyc.dataprep}
}
\details{
there probably isn't any reason to use this function on its own, but it can be optimized using "optim"}
\value{
returns the log likelihood of the model given data and parameter values}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Noah M. Reid}
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
