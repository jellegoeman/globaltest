\name{gt.result-class}
\docType{class}
\alias{result}
\alias{p.value}
\alias{gt.result-class}
\alias{show,gt.result-method}
\alias{result,gt.result-method}
\alias{p.value,gt.result-method}

\title{Class "gt.result" for results of the function globaltest}

\description{The class gt.result is the output of a call to
\code{\link{globaltest}} and the input of various plotting 
functions to visualize the test result. }

\section{Slots}{
  \describe{
    \item{\code{res}:}{Object of class "matrix". Test results, with a row for each element of 
    \code{test.genes}, and the following columns: \code{path.n}: length of input vector, 
    \code{test.n}: number of rows in \code{X} that matched the input vector, \code{Q}: test statistic, 
    \code{EQ}: expectation for the test statistic, \code{seQ}: standard deviation of the test 
    statistic, \code{p.val}: p value \code{comp.p}: comparative p if \code{\link[globaltest]{sampling = TRUE}}
    }
    \item{\code{X}:}{Object of class "matrix". The transformed data matrix.}
    \item{\code{Y}:}{Object of class "numeric". The transformed clinical outcome vector.}
    \item{\code{test.genes}:}{A list of vectors indicating the tested genes.}
    \item{\code{adjustmatrix}:}{Object of class "matrix" needed to calculate the expectation 
    of the test statistic.}
    \item{\code{Rsquare}:}{Object of class "numeric". Percentage of the variance of Y left after adjustment 
    for confounders or known risk factors.}
    \item{\code{model}:}{Object of class "character". The model used by the test ("linear" or "logistic").}
    \item{\code{df.adjust}:}{Object of class "numeric". The number of degrees of freedom used in the 
    adjustment of the test.} 
    }
}

\section{Methods}{ 
  \describe{
    \item{show}{(gt.result): Summarizes the test result}
    \item{p.value}{(gt.result): Extracts the p-values.}
    \item{result}{(gt.result): Extracts the results matrix.} 
    \item{geneplot}{(gt.result): Produces a plot to show the influence of individual
    genes on the test result produced by \code{\link{globaltest}}.}
    \item{sampleplot}{(gt.result): Produces a plot to show the influence of individual
    samples on the test result produced by \code{\link{globaltest}}.}
    \item{permutations}{(gt.result): Produces a histogram for visualization of the
    permutations used in \code{\link{globaltest}}.}
    \item{checkerboard}{(gt.result): Produces a plot to visualize the test result
    produced by \code{\link{globaltest}} by showing the association
    between pairs of samples.}
    \item{regressionplot}{(gt.result): Produces a plot which can be used to visualize
    the effect of specific samples on the test result produced by
    \code{\link{globaltest}}.}
    }
}

\author{Jelle Goeman: \email{j.j.goeman@lumc.nl}; Jan Oosting}

\seealso{\code{\link{globaltest}}, \code{\link{sampleplot}}, 
\code{\link{geneplot}}, \code{\link{permutations}}.}

\keyword{methods}
