\name{permutations}

\alias{permutations}

\title{Histogram for Global Test result}

\description{Produces a histogram for visualization of the permutations
used in \code{\link{globaltest}}.}

\usage{permutations(testresult)}

\arguments{
    \item{testresult}{The output of a call to \code{\link{globaltest}}.}}

\details{Produces a histogram of the value of the test statistic
calculated for permutations of the clinical outcome, compared with
the test statistic for the true clinical outcome.}

\note{Only works if the option \code{permutation = TRUE} was
chosen when calling \code{\link{globaltest}}.}

\value{\code{NULL} (no output).}

\references{J. J. Goeman, S.A. van de Geer, F. de Kort and J. C.
van Houwelingen, \emph{A global test for association of a group of
genes with a clinical outcome}, Technical Report MI 2003-03,
Mathematical Institute, Leiden University. Available from 
\url{http://www.math.leidenuniv.nl/~jgoeman}.}

\author{Jelle Goeman: \email{j.j.goeman@lumc.nl}}

\seealso{\code{\link{globaltest}}, \code{\link{checkerboard}},
\code{\link{geneplot}}, \code{\link{regressionplot}}.}

\examples{
    data(exampleX)      # Expression data (40 samples; 1000 genes)
    data(exampleY)      # Clinical outcome for the 40 samples
    pathway <- 1:25     # A pathway contains genes 1 to 25
    gt <- globaltest(exampleX, exampleY, test.genes = pathway, 
                                                permutation = TRUE)
    gt
    permutations(gt)
}

\keyword{hplot}
