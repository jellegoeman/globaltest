\name{gt.result-class}

\docType{class}

\alias{gt.result}
\alias(gt.result-class}

\title{Class gt.result, for results of the function globaltest}

\description{The class gt.result is the output of a call to
\code{\link{globaltest}} and the input of various plotting functions to
visualize the test reult.}

\section{Slots}{
\describe{
    \item{Q}{The test statistic.}
    \item{EQ}{The expected value of \code{Q} under the null
    hypothesis.}
    \item{seQ}{The standard deviation of \code{Q} under the null
    hypothesis.}
    \item{p.value}{The p-value resulting from the test.}
    \item{Qs}{Values of \code{Q} for permutations of the clinical
    variable}
    \item{matrixR}{The covariance matrix between the expression
    patterns of the samples.}
    \item{Y}{The centered clinical variable.}
    \item{influence}{The influence of each gene on the test
    statistic.}
    \item{exp.influence}{The expected influence of each gene on
    the test statistic under the null hypothesis.}
    \item{test.genes}{A logical vector indicating which genes were
    to be tested.}
    }
}

\section{Methods}{
\describe{
    \item{show}{Summarizes the test result}
    \item{checkerboard}{Produces a plot to visualize the test result
    produced by \code{\link{globaltest}} by showing the association
    between pairs of samples.}
    \item{geneplot}{Produces a plot to show the influence of individual
    genes on the test result produced by \code{\link{globaltest}}.}
    \item{permutations}{Produces a histogram for visualization of the
    permutations used in \code{\link{globaltest}}.}
    \item{regressionplot}{Produces a plot which can be used to visualize
    the effect of specific samples on the test result produced by
    \code{\link{globaltest}}.}
    }
}

\author{Jelle Goeman: \email{j.j.goeman@lumc.nl}}

\seealso{\code{\link{globaltest}}, \code{\link{regressionplot}},
\code{\link{geneplot}}, \code{\link{permutations}}, 
\code{\link{checkerboard}}.}

\keyword{methods}
