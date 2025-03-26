\name{checkrates}
\alias{checkrates}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
bGMYC output processing}
\description{
processes output from bgmyc.singlephy or bgmyc.multiphylo.}
\usage{
checkrates(result)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{result}{
output from bgmyc.singlephy or bgmyc.multiphylo}
}
\details{
this function takes the results of a bgmyc.singlephy or bgmyc.multiphylo run and outputs a matrix of all parameters, including the branching rates for the Yule and coalescent processes (this is the only place to get those last two). The branching rates are may be useful in evaluating the reliability of the model output (although I haven't tested this extensively). The GMYC is based on the idea that coalescent and speciation events happen at drastically different (perhaps on around an order of magnitude different) rates and thus can be differentiated in a phylogenetic tree. If the rates are similar, then the threshold point in the model might be relatively indistinct and the results unreliable. Comparison of the rates might help identify this. }
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{branchingrates }{a matrix containing bGMYC parameter estimates, including branching rates. lambda.div.mod is the effective Yule process branching parameter at the transition point between Yule and coalescent processes and thus of more interest than lambda.div.}
%% ...
}
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
