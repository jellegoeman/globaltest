\name{regressionplot}

\alias{regressionplot}

\title{Regression plot for Global Test result}

\description{Produces a plot which can be used to visualize the
effect of specific samples on the test result produced by
\code{\link{globaltest}}.}

\usage{regressionplot(testresult, samplenr = NULL)}

\arguments{
    \item{testresult}{The output of a call to \code{\link{globaltest}}.}
    \item{samplenr}{A vector giving row numbers of samples of interest.}}

\details{The regressionplot plots, for all pairs of samples, the
covariance between the expression patterns against the covariance
between their clinical outcomes. Each point in the plot therefore
represents a pair of samples. A regression line is fitted through
the samples, which visualizes the test result of the function
\code{\link{globaltest}}. A steeply increasing slope indicates a
high (possibly significant) value of the test statistic.

An optional argument \code{samplenr} can be supplied, giving
sample numbers of possibly outlying arrays. In this case, all
pairs of arrays involving one of the arrays in \code{samplenr} is
marked as a red cross, while the other pairs are marked as a blue
dot. The blue line which is fitted through all points can now be
compared to a red dotted line which is fitted though only the red
crosses.}

\value{\code{NULL} (no output).}

\references{J. J. Goeman, S.A. van de Geer, F. de Kort and J. C.
van Houwelingen, \emph{A global test for association of a group of
genes with a clinical outcome}, Technical Report MI 2003-03,
Mathematical Institute, Leiden University. Available from 
\url{http://www.math.leidenuniv.nl/~jgoeman}.}

\author{Jelle Goeman: \email{j.j.goeman@lumc.nl}}

\seealso{\code{\link{globaltest}}, \code{\link{checkerboard}},
  \code{\link{geneplot}}, \code{\link{permutations}}.} 

\examples{
if(interactive()){
    data(exampleX)      # Expression data (40 samples; 1000 genes)
    data(exampleY)      # Clinical outcome for the 40 samples
    pathway <- 1:25     # A pathway contains genes 1 to 25
    gt <- globaltest(exampleX, exampleY, test.genes = pathway)
    gt
    regressionplot(gt)
    regressionplot(gt,40)
}
}

\keyword{hplot}
