\name{checkerboard}

\alias{checkerboard}

\title{Checkerboard plot for Global Test result}

\description{Produces a plot to visualize the test result produced
by \code{\link{globaltest}}, by showing the association between
pairs of samples.}

\usage{checkerboard(gt, geneset = NULL, sort = TRUE, ...)}

\arguments{
    \item{gt}{The output of a call to \code{\link{globaltest}}.}
    \item{geneset}{The name or number of the geneset to be plotted 
(only necessary if multiple genesets were tested)}
    \item{sort}{A logical flag to indicate whether the samples should be sorted
by the clinical outcome to give a clearer picture.}
    \item{...}{Any extra arguments will be forwarded to the plotting function.}}

\details{The checkerboard shows the pairs of samples which have
high covariance in white and the pairs with low covariance in
black. This can be used to visualize the data and to search for
outlying arrays.}

\value{A matrix giving the old and the new sample numbers.}

\references{J. J. Goeman, S. A. van de Geer, F. de Kort and J. C. 
van Houwelingen, 2004, \emph{A global test for groups of genes: 
testing association with a clinical outcome}, 
\emph{Bioinformatics} 20 (1) 93--99. See also the How To 
Globaltest.pdf included with this package.}

\author{Jelle Goeman: \email{j.j.goeman@lumc.nl}; Jan Oosting}

\seealso{\code{\link{globaltest}}, \code{\link{regressionplot}},
\code{\link{geneplot}}, \code{\link{permutations}}.}

\examples{
if (interactive()) {
    data(exampleX)      # Expression data (40 samples; 1000 genes)
    data(exampleY)      # Clinical outcome for the 40 samples
    pathway <- 1:25     # A pathway contains genes 1 to 25
    gt <- globaltest(exampleX, exampleY, test.genes = pathway)
    gt
    checkerboard(gt)
}
}

\keyword{hplot}
