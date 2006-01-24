\name{geneplot}

\alias{geneplot}

\title{Gene Plot for Global Test}

\description{Produces a plot to show the influence of individual
genes on the test result produced by \code{\link{globaltest}}.}

\usage{geneplot(gt, geneset, genesubset, plot = TRUE, scale = FALSE,
            drawlabels = TRUE, labelsize = 0.6, ...)}

\arguments{
    \item{gt}{The output of a call to \code{\link{globaltest}}.}
    \item{geneset}{The name or number of the geneset to be plotted
(only necessary if multiple genesets were tested).}
    \item{genesubset}{A vector of names or numbers of genes to be plotted
(default: plot all genes)}
    \item{plot}{If \code{FALSE}: does not plot, but only returns a
\code{\link[gt.barplot-class]{gt.barplot}} object.}
    \item{scale}{Logical: should the bars be scaled to unit standard deviation?}
    \item{drawlabels}{Logical value to control drawing of gene ids (rownames) on the x-axis of the plot.}
    \item{labelsize}{Relative size of the labels on the x-axis. If it is \code{NULL} ,
the current value for \code{\link[graphics:par]{par("cex.axis")}}
is used}
    \item{...}{Any extra arguments will be forwarded to the plotting
function.}}

\details{The geneplot shows a bar an a reference line for each
gene. The bar shows the influence of each gene on the test
result: the test statistic Q is the average of the bars of the
genes. The reference line shows the expected influence if the
gene was not associated with the outcome. The marks on the bars
show the standard deviation of the bar under the null hypothesis.
The color of the bar indicates positive or negative correlation
of the gene with the clinical outcome, to distinguish between up
and downregulation.

The bottom margin is adjusted to allow enough space for the longest
gene id to draw under the axis.}

\value{An object of type \code{\link[gt.barplot-class]{gt.barplot}}.}

\references{For references, type: citation("globaltest"). See also
the vignette GlobalTest.pdf included with this package.}

\author{Jelle Goeman: \email{j.j.goeman@lumc.nl}; Jan Oosting}

\seealso{\code{\link{globaltest}}, \code{\link{sampleplot}}.}

\examples{
    # Breast cancer data (exprSet) from the Netherlands Cancer
    # Institute with annotation:
    data(vandeVijver)
    data(RosettaAnnotation)
    somePathways <- RosettaAnnotation[1:5]

    gt <- globaltest(vandeVijver, "StGallen", somePathways)

    if (interactive()){
      geneplot(gt[1])
    }

    gp <- geneplot(gt[1], plot = FALSE)
    if (interactive()){
      plot(sort(gp))
    }

}
\keyword{hplot}
