\name{sampling}

\alias{sampling}

\title{Sampling random 'pathways' for the Global Test}

\description{(DEPRECATED: please use use \code{\link{comparative}} instead). For every pathway in a result of \code{globaltest},
calculates how many randomly drawn groups of genes of the same size have a
smaller or equal p-value.}

\usage{sampling(gt, geneset, ndraws = 10^3)}

\arguments{
    \item{gt}{The output of a call to \code{\link{globaltest}}.}
    \item{geneset}{The name or number of the geneset(s) to be used
    (only necessary if multiple genesets were tested).}
    \item{ndraws}{The number of random pathways to be used.} }

\details{For every pathway in \code{gt[geneset]}, a number \code{ndraws}
random 'pathways' is selected by randomly sampling sets of genes of the same
size as the tested pathway. A 'comparative p-value' is calculated by counting
what proportion of the random pathways has a larger standardized test statistic
(Q - EQ) / sd(Q) than the tested pathway.}

\value{An object of class
\code{\link{gt.result}}.}

\note{
\code{sampling} has been deprecated. Please use \code{\link{comparative}} instead.

The function \code{sampling} cannot be applied to a
\code{\link{gt.result}} object resulting from a call to
\code{\link{permutations}}.}

\references{Jelle J. Goeman, Jan Oosting, Anne-Marie Cleton-Jansen, Jakob, K.
  Anninga, Hans C. van Houwelingen (2005) Testing association of a
  pathway with survival. Bioinformatics 21, 1950-1957. See also the vignette
Globaltest.pdf that comes with this package.}

\author{Jelle Goeman: \email{j.j.goeman@lumc.nl}; Jan Oosting}

\seealso{\code{\link{globaltest}}, \code{\link{permutations}},
\code{\link{sampleplot}}, \code{\link{geneplot}}.}

\examples{
    \dontshow{gt.options(warn.deprecated=FALSE)}
    # Breast cancer data (ExpressionSet) from the Netherlands Cancer
    # Institute with annotation:
    data(vandeVijver)
    data(annotation.vandeVijver)

    gt <- globaltest(vandeVijver, "StGallen", annotation.vandeVijver)

    sampling(gt, ndraws=100)
}

\keyword{htest}
