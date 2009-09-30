\name{gt.multtest}

\alias{gt.multtest}

\title{Correct globaltest results for multiple testing (Deprecated)}

\description{(DEPRECATED: please use use \code{\link{p.adjust}} instead). Corrects the raw p-values resulting from a call to
\code{globaltest} for multiple testing, using either Benjamini and
Hochberg's False Discovery Rate or Holm's procedure for controlling
the Family-Wise Error Rate.}

\usage{gt.multtest(gt, proc = c("FDR", "FWER"))}

\arguments{
    \item{gt}{The output of a call to \code{\link{globaltest}}.}
    \item{proc}{The procedure to be used. Either "FDR" for
    Benjamini and Hochberg's (1995) False Discovery Rate-controlling procedure or
    "FWER" for Holm's (1979) Family-Wise Error Rate controlling procedure.}
}

\details{This function is completely based on the
\code{\link[multtest:mt.rawp2adjp]{mt.rawp2adjp}} function from the
\code{multtest} package.}

\value{An object of class \code{\link{gt.result}}.}

\note{\code{gt.multtest} has been deprecated. Please use \code{\link{p.adjust}} instead.

This function must be called \emph{prior} to any selection of
significant genes.}

\references{For references, type: citation("globaltest"). See also
the vignette GlobalTest.pdf included with this package.}

\author{Jelle Goeman: \email{j.j.goeman@lumc.nl}; Jan Oosting}

\seealso{\code{\link{globaltest}}, \code{\link{gtGO}}.}

\examples{
    \dontshow{gt.options(warn.deprecated=FALSE)}
    # Breast cancer data (ExpressionSet) from the Netherlands Cancer
    # Institute with annotation:
    data(vandeVijver)
    data(annotation.vandeVijver)

    gt <- globaltest(vandeVijver, "StGallen", annotation.vandeVijver)
    sort(gt.multtest(gt))
}

\keyword{htest}
