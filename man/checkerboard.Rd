\name{checkerboard}

\alias{checkerboard}

\title{Checkerboard plot for Global Test result}

\description{Produces a plot to visualize the test result produced
by \code{\link{globaltest}}, by showing the association between
pairs of samples.}

\usage{checkerboard(testresult, sort = TRUE)}

\arguments{
    \item{testresult}{The output of a call to \code{\link{globaltest}}.}
    \item{sort}{A logical flag to indicate whether the samples should be sorted
by the clinical outcome to give a clearer picture.}}

\details{The checkerboard shows the pairs of samples which have
high covariance in white and the pairs with low covariance in
black. This can be used to visualize the data and to search for
outlying arrays.}

\value{A matrix giving the old and the new sample numbers.}

\references{J. J. Goeman, S.A. van de Geer, F. de Kort and J. C.
van Houwelingen, \emph{A global test for association of a group of
genes with a clinical outcome}, Technical Report MI 2003-03,
Mathematical Institute, Leiden University. Available from 
\url{http://www.math.leidenuniv.nl/~jgoeman}.}

\author{Jelle Goeman: \email{j.j.goeman@lumc.nl}}

\seealso{\code{\link{globaltest}}, \code{\link{regressionplot}},
\code{\link{geneplot}}, \code{\link{permutations}}.}

\examples{
    data(exampleX)      # Expression data (40 samples; 1000 genes)
    data(exampleY)      # Clinical outcome for the 40 samples
    pathway <- 1:25     # A pathway contains genes 1 to 25
    gt <- globaltest(exampleX, exampleY, test.genes = pathway)
    gt
    checkerboard(gt)
}

\keyword{hplot}
