\name{sampleplot}

\alias{sampleplot}

\title{Sample Plot for Global Test}

\description{Produces a plot to show the influence of individual
samples on the test result produced by \code{\link{globaltest}}.}

\usage{sampleplot(gt, geneset = NULL, samplesubset = NULL, drawlabels = TRUE, labelsize = 0.6, ...)}

\arguments{
    \item{gt}{The output of a call to \code{\link{globaltest}}.}
    \item{geneset}{The name or number of the geneset to be plotted 
(only necessary if multiple genesets were tested).}
    \item{samplesubset}{A vector of names or numbers of samples to be plotted 
(default: all samples)}
		\item{drawlabels}{Logical value to control drawing of the samplenames on the x-axis of the plot.}
		\item{labelsize}{Relative size of the labels on the x-axis. If it is \code{NULL} , 
		the current value for \code{\link[graphics:par]{par("cex.axis")}} is used}
    \item{...}{Any extra arguments will be forwarded to the plotting 
function.}}

\details{The sampleplot shows a bar and a reference line for each 
sample. The bar shows the influence of each gene on the test 
statistic. Samples with a positive influence carry evidence 
against the null hypothesis (in favour of a significant p-value), 
because they are are similar in expression profile to samples 
with a similar clinical outcome. Samples with a negative 
influence bar supply evidence in favour of the null hypothesis 
and of a non-significant p-value: they are relatively similar in 
expression profile to samples with a different clinical outcome.

The influence varies around zero if the tested geneset is not 
associated with the outcome. Marks on the bars show the 
standarddeviation of the influence under the null hypothesis for 
those samples which are more than one standard deviation away 
from zero.

The color of the bar indicates the sign of the residual of Y. In 
a logistic model the coloring this distinguishes the original 
groups.

The bottom margin is adjusted to allow enough space for the longest 
samplename to draw under the axis.}

\value{If samplenames were supplied in the call to 
\code{\link{globaltest}}, sampleplot returns the vector of 
samplenames corresponding to the samplenumbers appearing in the 
plot. Otherwise \code{NULL}.}

\references{J. J. Goeman, S. A. van de Geer, F. de Kort and J. C. 
van Houwelingen, 2004, \emph{A global test for groups of genes: 
testing association with a clinical outcome}, 
\emph{Bioinformatics} 20 (1) 93--99. See also the How To 
Globaltest.pdf included with this package.}

\author{Jelle Goeman: \email{j.j.goeman@lumc.nl}; Jan Oosting}

\seealso{\code{\link{globaltest}}, \code{\link{geneplot}}, 
\code{\link{regressionplot}}, \code{\link{checkerboard}}.}

\examples{
if (interactive()){ 
    data(exampleX)      # Expression data (40 samples; 1000 genes)
    data(exampleY)      # Clinical outcome for the 40 samples
    pathway <- 1:25     # A pathway contains genes 1 to 25
    gt <- globaltest(exampleX, exampleY, test.genes = pathway)
    gt
    sampleplot(gt)
}
}
\keyword{hplot}
