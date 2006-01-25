\name{globaltest}

\alias{globaltest}

\title{Global Test}

\description{In microarray data, tests a (list of) group(s) of
genes for significant association with a given clinical variable.}

\usage{globaltest(X, Y, genesets, model,
    levels, d, event = 1, adjust,
    method = c("auto", "asymptotic", "permutations", "gamma"),
    nperm = 10^4,
    scaleX = TRUE, ...) }

\arguments{
    \item{X}{Either a matrix of gene expression data, where columns correspond to
samples and rows to genes or a Bioconductor \code{\link[Biobase:exprSet-class]{exprSet}}. The data
should be properly normalized beforehand (and log- or otherwise
transformed), but missing values are allowed (coded as
\code{NA}). Gene and sample names can be included as the row and
column names of \code{X}.}
    \item{Y}{A vector with the clinical outcome of interest, having one value
for each sample. If \code{X} is an
\code{\link[Biobase:exprSet-class]{exprSet}} it can also be the
name of a covariate in the
\code{\link[Biobase:phenoData-class]{phenoData}} slot of the
exprSet, or a \code{\link[stats:formula]{formula}} object using
these names. If the clinical outcome is survival, \code{Y} should
contain the survival times.}
    \item{genesets}{Either a vector or a list of vectors. Indicates
the group(s) of genes to be tested. Each vector in
\code{genesets} can be given in three formats. Either it can be
a vector with 1 (\code{TRUE}) or 0 (\code{FALSE}) for each gene
in \code{X}, with 1 indicating that the gene belongs to the
group. Or it can be a vector containing the column numbers (in
\code{X}) of the genes belonging to the group. Or it can be a
subset of the rownames or \code{\link[Biobase:exprSet-class]{geneNames}}
for \code{X}.}
    \item{model}{Globaltest will try to determine the correct model from
the input of \code{Y} and \code{d}. To override the automatic choice, use
\code{model = "logistic"} for a two-valued outcome \code{Y} ,
\code{model = "linear"} for a continuous outcome and
\code{model = "survival"} for a survival outcome.}
    \item{levels}{If \code{Y} is a factor (or a category in the PhenoData slot of \code{X})
and contains more than 2 levels: \code{levels} is a vector of levels of \code{Y} to test. If
\code{levels} is length 2: test these 2 groups against each other.
If levels is length 1: test that level against the others.}
    \item{d}{A vector or the name of a covariate in the phenoData slot
of the exprSet \code{X}, to indicate which samples experienced an event.
Providing a value for \code{d} automatically sets \code{model = "survival"}}
    \item{event}{The value or values of \code{d} that indicates that
there was an event.}
    \item{adjust}{Confounders or risk factors for which the test must
be adjusted. Must be either a data frame or (if \code{X} is an
exprSet) the names of covariates in the phenoData slot \code{X}
or a \code{\link[stats:formula]{formula}} object using these
names. Default: no adjustment.}
    \item{method}{The method for calculation the p-value. Use code{method =
"asymptotic"} for the full asymptotic distribution of the test
statistic; \code{method = "gamma"} for the gamma (= scaled
chi-squared) approximation to that distribution and \code{method =
"permutations"} for a permutation p-value. The default: \code{method
= "auto"} chooses the permutations method if the number of possible
permutations does not exceed 10,000 and the asymptotic otherwise.
Note that \code{method = "gamma"} was the default option prior to
version 4.0.0.}
    \item{nperm}{A number of permutations. This gives the (maximum) number of
permutations to be used if \code{method = "permutations"} or
\code{method = "auto"}.}
    \item{scaleX}{If true, rescales the expression matrix to get pleasant value for
all test statistics. The expression matrix \code{X} is multiplied by a constant
in such a way that the expected value EQ of the test statistic for the global test
becomes exactly 10. This rescaling has no effect on the p-values.}
    \item{...}{Captures deprecated input for compatibility with older versions of
globaltest.} }

\details{The Global Test tests whether a group of genes (of any
size from one single gene to all genes on the array) is
significantly associated with a clinical variable. The group could
be for example a known pathway, an area on the genome or the set
of all genes. The test investigates whether samples with similar
clinical outcomes tend to have similar gene expression patterns.
For a significant result it is not necessary that the genes in
the group have similar expression patterns, only that many of
them are correlated with the outcome.}

\note{
The options globaltest options sampling and permutation
have been replaced by separate functions from version 3.0. See
\code{\link{sampling}} and \code{\link{permutations}}.}

\value{The function returns an object of class
\code{\link[gt.result-class]{gt.result}}.}

\references{For references, type: \code{citation("globaltest")}. See
also the vignette GlobalTest.pdf included with this package.}

\author{Jelle Goeman: \email{j.j.goeman@lumc.nl}; Jan Oosting}

\seealso{Many more examples in the vignette! \code{\link{geneplot}},
\code{\link{sampleplot}}, \code{\link{sampling}},
\code{\link{gt.multtest}}, \code{\link{permutations}},
\code{\link{checkerboard}}, \code{\link{regressionplot}}.}

\examples{
    # Breast cancer data (exprSet) from the Netherlands Cancer
    # Institute with annotation:
    data(vandeVijver)
    data(annotation.vandeVijver)

    # Many possible calls. See the vignette for more examples and explanation.
    globaltest(vandeVijver, "StGallen")
    globaltest(vandeVijver, "StGallen", annotation.vandeVijver)
    globaltest(vandeVijver, "Surv(TIMEsurvival, EVENTdeath)", annotation.vandeVijver)
    globaltest(vandeVijver, StGallen ~ Posnodes + StGallen, annotation.vandeVijver)
    globaltest(vandeVijver, "StGallen", method = "p")

    # Store the test result
    # See help(gt.result) for more options
    gt <- globaltest(vandeVijver, "StGallen", annotation.vandeVijver)
    gt[1:2]
    sort(gt)
    p.value(gt)

    # Also with simple vector/matrix input
    X <- matrix(rnorm(3000), 100, 30)  # random expression data
    Y <- 1:30                          # a response variable
    pathway <- 1:40                    # a pathway

    globaltest(X, Y)
    globaltest(X, Y, pathway)
}

\keyword{htest}
