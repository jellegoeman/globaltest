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

\details{A permutation p-value is calculated by comparing
the value of the test statistic for random permutations of the clinical outcome
to the value of the test statistic for the true clinical outcome.
The permuted test statistics are stored for later visualisation using
\code{hist}.}

\value{An object of class
\code{\link[gt.result-class]{gt.result}}.}

\note{Permutations does not work if the adjusted version of
globaltest was used.}

\references{J. J. Goeman, S. A. van de Geer, F. de Kort and J. C.
van Houwelingen, 2004, \emph{A global test for groups of genes:
testing association with a clinical outcome},
\emph{Bioinformatics} 20 (1) 93--99. See also the How To
Globaltest.pdf included with this package.}

\author{Jelle Goeman: \email{j.j.goeman@lumc.nl}; Jan Oosting}

\seealso{\code{\link{globaltest}}, \code{\link{sampleplot}},
\code{\link{geneplot}}.}

\examples{
    data(exampleX)      # Expression data (40 samples; 1000 genes)
    data(exampleY)      # Clinical outcome for the 40 samples
    pathway <- 1:25     # A pathway contains genes 1 to 25

    gt <- globaltest(exampleX, exampleY, test.genes = pathway)
    gt

    perm.gt <- permutations(gt)
    perm.gt

    if (interactive()) {
      hist(perm.gt)
    }
}
\keyword{htest}
