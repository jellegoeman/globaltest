\name{permutations}

\alias{permutations}

\title{Permutations Plot for Global Test}

\description{Produces a histogram for visualization of the permutations
used in \code{\link{globaltest}}.}

\usage{permutations(gt, geneset = NULL, nperm = 10^4)}

\arguments{
    \item{gt}{The output of a call to \code{\link{globaltest}}.}
    \item{geneset}{The name or number of the geneset to be plotted 
(only necessary if multiple genesets were tested).}
    \item{nperm}{The number of permutations to be used.} }

\details{Produces a histogram of the value of the test statistic
calculated for permutations of the clinical outcome, compared with
the test statistic for the true clinical outcome.}

\value{\code{NULL} (no output except the plot).}

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
if (interactive()) {
    data(exampleX)      # Expression data (40 samples; 1000 genes)
    data(exampleY)      # Clinical outcome for the 40 samples
    pathway <- 1:25     # A pathway contains genes 1 to 25
    gt <- globaltest(exampleX, exampleY, test.genes = pathway, 
                                                permutation = TRUE)
    gt
    permutations(gt)
}
}
\keyword{hplot}
