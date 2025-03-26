\name{bgmyc.point}
\alias{bgmyc.point}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
returns a point estimate of species limits based on a posterior probability threshold input by the user}
\description{
this function processes output from a bgmyc analysis and produces a single point estimate of species limits based on a posterior probability threshold}
\usage{
bgmyc.point(probmat, ppcutoff)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{probmat}{
output from function spec.probmat of class "bgmycprobmat"}
  \item{ppcutoff}{
a posterior probability threshold for lumping samples into species. e.g. if 0.05 were selected, all individuals having greater than 0.05 posterior probability of conspecificity will be lumped into species returned.}
}

\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
returns an R list where each entry gives the members of a species that meets the threshold. it is worth noting that sometimes, because of uncertainty in the phylogenetic tree, not all individuals in a "species" resulting from a given threshold will have greater than the specified probability of conspecificity. they can all be connected, however, through individuals whose probability of conspecificity is greater than the threshold.}
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
