\name{geneplot}

\alias{geneplot}

\title{Geneplot for Global Test result}

\description{Produces a plot to show the influence of individual
genes on the test result produced by \code{\link{globaltest}}.}

\usage{geneplot(testresult)}

\arguments{
    \item{testresult}{The output of a call to \code{\link{globaltest}}.}}

\details{The geneplot shows a bar an a reference line for each
gene. The bar shows the influence of each gene on the test result;
the reference line shows the expected influence if the gene was
not associated with the outcome. The color of the bar indicates
positive or negative correlation of the gene with the clinical
outcome, to distinguish between up and downregulation.}

\value{If gene names were supplied in the call to 
\code{\link{globaltest}}, geneplot returns a data.frame referring 
the gene names to the gene numbers appearing in the plot. 
Otherwise \code{NULL}.}

\references{J. J. Goeman, S.A. van de Geer, F. de Kort and J. C. 
van Houwelingen, \emph{A global test for association of a group 
of genes with a clinical outcome}, Technical Report MI 2003-03, 
Mathematical Institute, Leiden University. Available from 
\url{http://www.math.leidenuniv.nl/~jgoeman}.}

\author{Jelle Goeman: \email{j.j.goeman@lumc.nl}}

\seealso{\code{\link{globaltest}}, \code{\link{regressionplot}},
  \code{\link{checkerboard}}, \code{\link{permutations}}.}  

\examples{
if(interactive()){
    data(exampleX)      # Expression data (40 samples; 1000 genes)
    data(exampleY)      # Clinical outcome for the 40 samples
    pathway <- 1:25     # A pathway contains genes 1 to 25
    gt <- globaltest(exampleX, exampleY, test.genes = pathway)
    gt
        geneplot(gt)
}
}

\keyword{hplot}
