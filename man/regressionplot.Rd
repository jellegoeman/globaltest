\name{regressionplot}

\alias{regressionplot}

\title{Regression Plot for Global Test}

\description{Produces a plot which can be used to visualize the
effect of specific samples on the test result produced by
\code{\link{globaltest}}.}

\usage{regressionplot(gt, geneset = NULL, sampleid = NULL, ...)}

\arguments{
    \item{gt}{The output of a call to \code{\link{globaltest}}.}
    \item{geneset}{The name or number of the geneset to be plotted 
(only necessary if multiple genesets were tested).}
    \item{sampleid}{A vector of names or numbers of the samples of interest.}
    \item{...}{Any extra arguments will be forwarded to the plotting function.}}

\details{The regressionplot plots, for all pairs of samples, the
covariance between the expression patterns against the covariance
between their clinical outcomes. Each point in the plot therefore
represents a pair of samples. A regression line is fitted through
the samples, which visualizes the test result of the function
\code{\link{globaltest}}. A steeply increasing slope indicates a
high (possibly significant) value of the test statistic.

An optional argument \code{sampleid} can be supplied, giving
sample numbers of possibly outlying arrays. In this case, all
pairs of arrays involving one of the arrays in \code{sampleid} is
marked as a red cross, while the other pairs are marked as a blue
dot. The blue line which is fitted through all points can now be
compared to a red dotted line which is fitted though only the red
crosses.}

\value{\code{NULL} (no output).}

\note{Regressionplot does not work if the adjusted version of 
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
    gt <- globaltest(exampleX, exampleY, test.genes = pathway)
    gt
    regressionplot(gt)
    regressionplot(gt, sampleid = 40)
}
}
\keyword{hplot}
