\name{checkerboard}

\alias{checkerboard}

\title{Checkerboard plot for Global Test (Deprecated)}

\description{Produces a plot to visualize the test result produced
by \code{\link{globaltest}}, by showing the association between
pairs of samples.}

\usage{checkerboard(gt, geneset, sort = TRUE, drawlabels = TRUE,
labelsize = 0.6, ...)}

\arguments{
    \item{gt}{The output of a call to \code{\link{globaltest}}.}
    \item{geneset}{The name or number of the geneset to be plotted
(only necessary if multiple genesets were tested)}
    \item{sort}{A logical flag to indicate whether the samples should be sorted
by the clinical outcome to give a clearer picture.}
        \item{drawlabels}{Logical value to control drawing of the samplenames on the x- and y-axis of the plot.}
        \item{labelsize}{Relative size of the labels on the x- and y-axis. If it is \code{NULL} ,
        the current value for \code{\link[graphics:par]{par("cex.axis")}} is used}
    \item{...}{Any extra arguments will be forwarded to the plotting function.}}

\details{The checkerboard shows the pairs of samples which have
high covariance in white and the pairs with low covariance in
black. This can be used to visualize the data and to search for
outlying arrays.

The left and bottom margins are adjusted to allow enough space
for the longest samplename.}

\value{A matrix giving the old and the new sample numbers.}

\note{\code{checkerboard} has been deprecated. Please use \code{\link{subjects}} instead.}


\references{J. J. Goeman, S. A. van de Geer, F. de Kort and J. C.
van Houwelingen, 2004, \emph{A global test for groups of genes:
testing association with a clinical outcome},
\emph{Bioinformatics} 20 (1) 93--99. See also the How To
Globaltest.pdf included with this package.}

\author{Jelle Goeman: \email{j.j.goeman@lumc.nl}; Jan Oosting}

\seealso{\code{\link{globaltest}}, \code{\link{sampleplot}},
\code{\link{geneplot}}.}

\examples{
    # Breast cancer data (ExpressionSet) from the Netherlands Cancer
    # Institute with annotation:
    data(vandeVijver)
    data(annotation.vandeVijver)
    aPathway <- annotation.vandeVijver[1]

    gt <- globaltest(vandeVijver, "ESR1", aPathway)

    if (interactive()){
      checkerboard(gt[1])
    }
}

\keyword{hplot}
