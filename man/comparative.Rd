\name{Comparative proportions}

\alias{comparative}

\title{Comparative proportions for the Global Test}

\description{Comparing the result of a global test performed on a subset to the test results of random subsets of the same size.}

\usage{
comparative(object, N = 1e3, z.scores = TRUE, trace)
}

\arguments{
    \item{object}{A \code{\link{gt.object}}, usually one in which one or more subsets of a large number of covariates were tested.}
    \item{N}{The number of random subsets to generate.}
    \item{z.scores}{If set to \code{TRUE}, compares the subset to random subsets on the basis of the z-scores of the test. If \code{FALSE}, uses the p-values instead.}
    \item{trace}{If set to \code{TRUE}, reports progress information. The default is to set trace to \code{TRUE} if R is in \code{\link{interactive}} mode and more than one comparative proportion is to be calculated.}
}

\details{In a situation when many covariates out of a large set are associated with the response, it is sometimes interesting to know p-value of the tested subset compares to random subsets of the same size. The \code{comparative} function calculates the proportion of random subsets of the covariates of the same size as the tested subset that have a better score than the tested subset. This proportion is a diagnostic tool to help interpret the test result; it should not be interpreted as a p-value.
}

\value{An object of class \code{\link{gt.object}} with an appropriate column added to the test results matrix.}

\references{
The comparative proportion is an enrichment type analysis. For the pros and cons of such an analysis, see

Goeman and Buhlmann (2007) Analyzing gene expression data in terms of gene sets: methodological issues. Bioinformatics 23 (8) 980-987.
}

\author{Jelle Goeman: \email{j.j.goeman@lumc.nl}; Jan Oosting}

\seealso{
The \code{\link{gt}} function. The \code{\link{gt.object}} function and useful functions associated with that object.
}

\examples{
    # Simple examples with random data here
    # Real data examples in the Vignette

    # Random data: covariates A,B,C are correlated with Y
    set.seed(1)
    Y <- rnorm(20)
    X <- matrix(rnorm(200), 20, 10)
    X[,1:3] <- X[,1:3] + Y
    colnames(X) <- LETTERS[1:10]

    # Some subsets of interest
    my.sets <- list(c("A", "B"), c("C","D"), c("D", "E"))

    # Comparative proportions
    res <- gt(Y, X, subsets = my.sets)
    comparative(res)
    comparative(res, z.scores=FALSE)
}

\keyword{htest} 