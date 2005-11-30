\name{globaltest}

\alias{globaltest}

\title{Global Test}

\description{In microarray data, tests a (list of) group(s) of
genes for significant association with a given clinical variable.}

\usage{globaltest(X, Y, test.genes, model,
    levels, d, event = 1, adjust, scaleX = TRUE, ...) }

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
    \item{test.genes}{Either a vector or a list of vectors. Indicates
the group(s) of genes to be tested. Each vector in
\code{test.genes} can be given in three formats. Either it can be
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
1. The options globaltest options sampling and permutation
have been replaced by separate functions from version 3.0. See
\code{\link{sampling}} and \code{\link{permutations}}.

2. The scaling of the test statistic Q is arbitrary and does not
influence the test result. The scaling used in the software
slightly different from the unscaled version presented in the paper.
The mechanism of the rescaling is as follows: Remember that Q is a
average of a Qi for all genes i in the tested geneset, so EQ is the
average of the EQi. The rescaling is such that the EQ for the pathway
of all genes (i.e. the mean of all EQi) is 10.}

\value{The function returns an object of class
\code{\link[gt.result-class]{gt.result}}.}

\references{J. J. Goeman, S. A. van de Geer, F. de Kort and J. C.
van Houwelingen, 2004, \emph{A global test for groups of genes:
testing association with a clinical outcome},
\emph{Bioinformatics} 20 (1) 93--99. See also the vignette
Globaltest.pdf included with this package.}

\author{Jelle Goeman: \email{j.j.goeman@lumc.nl}; Jan Oosting}

\seealso{\code{\link{geneplot}}, \code{\link{sampleplot}},
\code{\link{sampling}}, \code{\link{permutations}},
\code{\link{checkerboard}}, \code{\link{regressionplot}}.}

\examples{
    data(exampleX)      # Expression data (40 samples; 1000 genes)
    data(exampleY)      # Clinical outcome for the 40 samples
    pathway1 <- 1:25    # A pathway contains genes 1 to 25
    pathway2 <- 26:50   # another pathway
    gt <- globaltest(exampleX, exampleY, list(pathway1,pathway2))
    gt
}
\keyword{htest}
