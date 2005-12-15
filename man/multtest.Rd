\name{multtest}

\alias{multtest}

\title{Correct globaltest results for multiple testing}

\description{Corrects the raw p-values resulting from a call to
\code{globaltest} for multiple testing, using either Benjamini and
Hochberg's False Discovery Rate or Holm's procedure for controlling
the Family-Wise Error Rate.}

\usage{multtest(gt, proc = c("FDR", "FWER"))}

\arguments{
    \item{gt}{The output of a call to \code{\link{globaltest}}.}
    \item{proc}{The procedure to be used. Either "FDR" for
    Benjamini and Hochberg's (1995) False Discovery Rate-controlling procedure or
    "FWER" for Holm's (1979) Family-Wise Error Rate controlling procedure.}
}

\details{This function is completely based on the
\code{\link[multtest:mt.rawp2adjp]{mt.rawp2adjp}} function from the
\code{multtest} package.}

\value{An object of class \code{\link[gt.result-class]{gt.result}}.}

\note{This function must be called \emph{prior} to any selection of
significant genes.}

\references{For references, type: citation("globaltest"). See also
the vignette GlobalTest.pdf included with this package.}

\author{Jelle Goeman: \email{j.j.goeman@lumc.nl}; Jan Oosting}

\seealso{\code{\link{globaltest}}.}

\examples{
    data(exampleX)      # Expression data (40 samples; 1000 genes)
    data(exampleY)      # Clinical outcome for the 40 samples
    pathway1 <- 1:25    # A pathway contains genes 1 to 25
    pathway2 <- 26:50   # another pathway
    gt <- globaltest(exampleX, exampleY, list(pathway1,pathway2))

    multtest(gt)
}

\keyword{htest}
