\name{gt.object class}
\docType{class}
\alias{p.value}
\alias{gt.object}
\alias{result}
\alias{size}
\alias{gt.object-class}
\alias{show,gt.object-method}
\alias{summary,gt.object-method}
\alias{size,gt.object-method}
\alias{[,gt.object-method}
\alias{[[,gt.object-method}
\alias{length,gt.object-method}
\alias{weights,gt.object-method}
\alias{hist,gt.object-method}
\alias{result,gt.object-method}
\alias{model.matrix,gt.object-method}
\alias{p.value,gt.object-method}
\alias{z.score,gt.object-method}
\alias{names,gt.object-method}
\alias{names<-,gt.object-method}
\alias{alias,gt.object-method}
\alias{alias<-}
\alias{alias<-,gt.object-method}
\alias{sort,gt.object-method}


\title{Class "gt.object" for storing the result of the function gt}

\description{The class gt.object is the output of a call to
\code{\link{gt}}. It stores the information needed for various diagnostic plots. }

\section{Slots}{
  These slots are not meant to be directly accessed by the user.
  \describe{
    \item{\code{result}:}{Object of class "matrix". The number of rows of this matrix is the number of tests performed. The matrix has at least the columns "p-value", "Statistic" "Expected", "Std.dev", and "#Cov".}
    \item{\code{extra}:}{Object of class "data.frame". Holds additional information that may be added later about the tests performed, such as multiplicity-adjusted p-values (see \code{\link{multtest}}), alias names for tests and comparative proportions (see \code{\link{comparative}}).}
    \item{\code{call}:}{The matched call to \code{\link{gt}}.}
    \item{\code{functions}:}{A "list" of various functions used by the \code{\link{covariates}} and \code{\link{subjects}} functions and the various methods.}
    \item{\code{subsets}:}{A "list" or "NULL". Stores the subsets tested, if more than one.}
    \item{\code{structure}:}{A "list" or "NULL". Stores subset and superset relationships between the sets in the "subsets" slot.}
    \item{\code{weights}:}{A "list" or "NULL". Stores the weight vectors used for testing, if more than one.}
    \item{\code{alternative}:}{If \code{\link{gt}} was called with \code{x = TRUE}, stores the design matrix of the alternative hypothesis; "NULL" otherwise.}
    \item{\code{null}:}{If \code{\link{gt}} was called with \code{x = TRUE}, stores the design matrix of the null hypothesis; "NULL" otherwise.}
    \item{\code{directional}}{Stores the \code{directional} argument of the call to \code{\link{gt}}.}
    \item{\code{legend}}{Object of class "list". Stores appropriate legends for the \code{\link{covariates}} and \code{\link{subjects}} plots.}
    \item{\code{model}}{Object of class "character". Stores the model.}
  }
}

\section{Methods}{
  \describe{
    \item{show}{(gt.object): Prints the test results: p-value, test statistic, expected value of the test statistic under the null hypothesis, standard deviation of the test statistic under the null hypothesis, and number of covariates tested.}
    \item{summary}{(gt.object): Prints the test results (as \code{show}) plus additional information on the model and the test.}
    \item{p.value}{(gt.object): Extracts the p-values.}
    \item{z.score}{(gt.object): Extracts z-score: (Test statistic - Expected value) / Standard deviation.}
    \item{result}{(gt.object): Extracts the results matrix together with the additional (e.g. multiple testing) information in the \code{multtest} slot.}
    \item{sort}{(gt.object): Sorts the pathways to increasing p-values. Equal p-values are sorted on decreasing z-scores.}
    \item{"["}{(gt.object): Extracts results of one or more test results if multiple tests were performed. Identical to "[[".}
    \item{"[["}{(gt.object): Extracts results of one or more test results if multiple tests were performed. Identical to "[".}
    \item{length}{(gt.object): The number of tests performed.}
    \item{size}{(gt.object): Extracts a vector with the number of alternative covariates tested for each test.}
    \item{names}{(gt.object): Extracts the row names of the results matrix.}
    \item{names<-}{(gt.object): Changes the row names of the results matrix. Duplicate names are not allowed, but see \code{alias}.}
    \item{alias}{(gt.object): Extracts the "alias" column of the results matrix that can be used to add additional information on each test perfomed.}
    \item{alias<-}{(gt.object): Changes the "alias" column of the results matrix. Note that unlike for names, duplicate aliases are allowed.}
    \item{weights}{(gt.object): extracts the effective weights of the covariates as they are used internally by the test.}
    \item{hist}{(gt.object): Produces a histogram to visualize the permutation test statistics. Only relevant after permutation testing.}
    \item{covariates}{(gt.object): Produces a plot to show the influence of individual covariates on the test result. See \code{\link{covariates}} for details.}
    \item{subjects}{(gt.object): Produces a plot to show the influence of individual subjects on the test result. See \code{\link{subjects}} for details.}
    \item{multtest}{(gt.object): Performs multiple testing correction and produces multiplicity-corrected p-values. See \code{\link{multtest}} for details.}
    \item{comparative}{(gt.object): Compares the p-values of tests performed on a subsets or weights with p-values of random subsets of covariates of same size or randomly distributed weights.  See \code{\link{comparative}} for details.}
    }
}

\author{Jelle Goeman: \email{j.j.goeman@lumc.nl}; Jan Oosting}

\seealso{\code{\link{gt}}, \code{\link{covariates}}, \code{\link{subjects}}.}

\keyword{methods}


\examples{
    # Simple examples with random data here
    # Real data examples in the Vignette

    # Random data: covariates A,B,C are correlated with Y
    Y <- rnorm(20)
    X <- matrix(rnorm(200), 20, 10)
    X[,1:3] <- X[,1:3] + 0.5*Y
    colnames(X) <- LETTERS[1:10]

    # Make a gt.object
    sets <- list(odd = c(1,3,5,7,9), even = c(2,4,6,8,10))
    res <- gt(Y, X, subsets=sets)

    # Show the results
    res
    summary(res)
    sort(res)
    p.value(res)

    # Names
    names(res)
    names(res) <- c("ODD", "EVEN")
    alias(res) <- c("odd covariates", "even covariates")

    # Multiple testing
    multtest(res, method = "Holm")
    multtest(res, method = "BH")

    # Diagnostics
    weights(res)
    covariates(res[1])
    subjects(res[1])

    # Permutation testing
    res <- gt(Y, X, perm = 1e4)
    hist(res)

} 