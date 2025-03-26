\name{bGMYC-package}
\alias{bGMYC-package}
\alias{bGMYC}
\docType{package}
\title{
A Bayesian MCMC implementation of the general mixed Yule-coalescent model for species delimitation. 
}
\description{This package implements a Bayesian version of the general mixed Yule-coalescent model for species delimitation introduced by Pons et al. 2006. It uses Markov Chain Monte Carlo simulation to evaluate the posterior distribution of the model and uses a flexible method of specifying prior distributions on the model parameters. The primary strength of this method over the pre-existing implementation is that it allows the use of multiple trees sampled from an MCMC analysis, eliminating the conditioning of results on a single point estimate of the phylogenetic tree, which often is associated with high uncertainty. 
}
\details{
\tabular{ll}{
Package: \tab bGMYC\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2012-07-26\cr
License: \tab GPLv3\cr
}
}
\author{
Noah M. Reid
Maintainer: Noah M. Reid <noah.reid@gmail.com>
}
\references{
Pons et al. 2006
Reid and Carstens 2012
}
\examples{
#for user data, ape tree format is required. easiest to input a nexus tree file using ape function "read.nexus"
#in this case, load test data, an ape multiPhylo object
data(east10)

#two primary functions for running analysis: bgmyc.singlephy and bgmyc.multiphylo
#run bGMYC on a single tree (phylo object)
result.single<-bgmyc.singlephy(phylo=east10[[1]], mcmc=100, burnin=1, thinning=10)

#run bGMYC on multiple trees (multiphylo object)
result.multi<-bgmyc.multiphylo(multiphylo=east10, mcmc=100, burnin=1, thinning=10)

\dontrun{
#several functions for processing and visualizing output
#visualize MCMC output
plot(result.single)
}

#for multiple trees, visualize MCMC output OR pool samples from across trees into a single matrix using plot=FALSE
plot(result.multi, plot=FALSE)->parameter.matrix

#get the posterior probabilities of each possible species cluster
#optionally writes to an output file. 
bgmyc.spec(res=result.multi)->result.spec

#produce a sequence by sequence matrix of probabilities of conspecificity
spec.probmat(result.multi)->result.probmat

\dontrun{
#plot the matrix of conspecificity probabilities
#provide a single tree to order the matrix
plot(result.probmat, east10[[1]])
}



}
