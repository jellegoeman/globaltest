\name{gt.result-class}
\docType{class}
\alias{gt.result-class}
\alias{show,gt.result-method}

\title{Class "gt.result" for results of the function globaltest}

\description{The class gt.result is the output of a call to
\code{\link{globaltest}} and the input of various plotting 
functions to visualize the test result.}

\section{Methods}{ \describe{
    \item{show}{Summarizes the test result}
    \item{p.value}{Extracts the p-values.}
    \item{result}{Extracts the results matrix.} 
    \item{geneplot}{Produces a plot to show the influence of individual
    genes on the test result produced by \code{\link{globaltest}}.}
    \item{sampleplot}{Produces a plot to show the influence of individual
    samples on the test result produced by \code{\link{globaltest}}.}
    \item{permutations}{Produces a histogram for visualization of the
    permutations used in \code{\link{globaltest}}.}
    \item{checkerboard}{Produces a plot to visualize the test result
    produced by \code{\link{globaltest}} by showing the association
    between pairs of samples.}
    \item{regressionplot}{Produces a plot which can be used to visualize
    the effect of specific samples on the test result produced by
    \code{\link{globaltest}}.}
    }
}

\section{Slots}{
\describe{
    \item{res}{The matrix of test results.}
    \item{X}{The transformed data matrix.}
    \item{Y}{The transformed outcome vector.}
    \item{test.genes}{A list of vectors indicating the tested genes.}
    \item{adjustmatrix}{A matrix needed to calculate the expectation 
    of the test statistic.}
    \item{Rsquare}{Percentage of the variance of Y left after adjustment 
    for confounders or known risk factors.}
    \item{model}{The model used by the test (linear or logistic).}
    \item{df.adjust}{The number of degrees of freedom used in the 
    adjustment of the test.} 
    }
}

\author{Jelle Goeman: \email{j.j.goeman@lumc.nl}; Jan Oosting}

\seealso{\code{\link{globaltest}}, \code{\link{p.value}}, 
\code{\link{result}}, \code{\link{sampleplot}}, 
\code{\link{geneplot}}, \code{\link{permutations}}.}

\keyword{classes}
