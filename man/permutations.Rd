\name{permutations}

\alias{permutations}

\title{Permutations version of the Global Test}

\description{Produces permutations of the \code{\link{globaltest}} test
statistic and uses these to recalculate the p-value.}

\usage{permutations(gt, geneset, nperm = 10^4)}

\arguments{
    \item{gt}{The output of a call to \code{\link{globaltest}}.}
    \item{geneset}{The name or number of the geneset(s) to be used
    (only necessary if multiple genesets were tested).}
    \item{nperm}{The number of permutations to be used.} }

\details{A permutation p-value is calculated by comparing the value
of the test statistic using permutations of the clinical outcome to
the value of the test statistic for the true clinical outcome. If
the number of possible permutations is smaller than the value of
\code{nperm}, all possible permutations are used, otherwise
\code{nperm} random permutations. The permuted test statistics are
stored for later visualisation using \code{\link{hist}}.}

\value{An object of class
\code{\link[gt.result-class]{gt.result}}.}

\note{Permutations does not work if the adjusted version of
globaltest was used.}

\references{For references, type: citation("globaltest"). See also
the vignette GlobalTest.pdf included with this package.}

\author{Jelle Goeman: \email{j.j.goeman@lumc.nl}; Jan Oosting}

\seealso{\code{\link{globaltest}}, \code{\link{sampleplot}},
\code{\link{geneplot}}.}

\examples{
    # Breast cancer data (exprSet) from the Netherlands Cancer
    # Institute with annotation:
    data(vandeVijver)
    data(annotation.vandeVijver)
    aPathway <- annotation.vandeVijver[1]

    gt <- globaltest(vandeVijver, "StGallen", aPathway)
    perm.gt <- permutations(gt)

    if (interactive()) {
      hist(perm.gt)
    }
}
\keyword{htest}
