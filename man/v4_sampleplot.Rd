\name{sampleplot}

\alias{sampleplot}

\title{Sample Plot for Global Test (Deprecated)}

\description{Produces a plot to show the influence of individual
samples on the test result produced by \code{\link{globaltest}}.}

\usage{sampleplot(gt, geneset, samplesubset, scale = TRUE, drawlabels = TRUE,
  labelsize = 0.6, plot = TRUE, addlegend = TRUE, ...)}

\arguments{
    \item{gt}{The output of a call to \code{\link{globaltest}}.}
    \item{geneset}{The name or number of the geneset to be plotted
      (only necessary if multiple genesets were tested).}
    \item{samplesubset}{A vector of names or numbers of samples to be plotted.
      Default: plot all samples}
    \item{scale}{Logical: should the bars be scaled to unit standard deviation?}
    \item{drawlabels}{Logical value to control drawing of the samplenames on the x-axis of the plot.}
    \item{labelsize}{Relative size of the labels on the x-axis. If it is \code{NULL} ,
        the current value for \code{\link[graphics:par]{par("cex.axis")}} is used}
    \item{plot}{If \code{FALSE}: does not plot, but only returns a
      \code{\link[gt.barplot-class]{gt.barplot}} object.}
    \item{addlegend}{If \code{FALSE}: does not add a legend to the plot
        or to the \code{\link[gt.barplot-class]{gt.barplot}} object.}
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

\value{An object of type \code{\link[gt.barplot-class]{gt.barplot}}.}

\note{\code{sampleplot} has been deprecated. Please use \code{\link{subjects}} instead.}

\author{Jelle Goeman: \email{j.j.goeman@lumc.nl}; Jan Oosting}

\seealso{\code{\link{globaltest}}, \code{\link{geneplot}},
\code{\link{regressionplot}}, \code{\link{checkerboard}}.}

\examples{
    \dontshow{gt.options(warn.deprecated=FALSE)}
    # Breast cancer data (ExpressionSet) from the Netherlands Cancer
    # Institute with annotation:
    data(vandeVijver)
    data(annotation.vandeVijver)

    gt <- globaltest(vandeVijver, "StGallen", annotation.vandeVijver)

    if (interactive()){
      sampleplot(gt[1])
    }

    sp <- sampleplot(gt[1], plot = FALSE)
    if (interactive()){
      plot(sort(sp))
    }
}
\keyword{hplot}
