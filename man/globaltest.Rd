\name{globaltest}

\alias{globaltest}

\title{Global Test}

\description{In microarray data, tests a (list of) group(s) of 
genes for significant association with a given clinical variable.}

\usage{globaltest(X, Y, test.genes = NULL, model = NULL, 
    levels = NULL, adjust = NULL, permutation = FALSE, nperm = NULL, 
    sampling = FALSE, ndraws = NULL, verbose = TRUE) }

\arguments{
    \item{X}{Either a matrix of gene expression data, where columns correspond to
samples and rows to genes or a Bioconductor \code{\link[Biobase:exprSet-class]{exprSet}}. The data 
should be properly normalized beforehand (and log- or otherwise 
transformed), but missing values are allowed (coded as 
\code{NA}). Gene and sample names can be included as the row and 
column names of \code{X}.}
    \item{Y}{A vector with the clinical outcome of interest, having one value
for each sample. If X is an \code{\link[Biobase:exprSet-class]{exprSet}} it can 
also be the name of a covariate in the \code{\link[Biobase:phenoData-class]{phenoData}} slot of the exprSet}
    \item{test.genes}{Either a vector or a list of vectors. Indicates 
the group(s) of genes to be tested. Each vector in 
\code{test.genes} can be given in three formats. Either it can be 
a vector with 1 (\code{TRUE}) or 0 (\code{FALSE}) for each gene 
in \code{X}, with 1 indicating that the gene belongs to the 
group. Or it can be a vector containing the column numbers (in 
\code{X}) of the genes belonging to the group. Or it can be a 
subset of the rownames or \code{\link[Biobase:exprSet-class]{geneNames}} for \code{X}.}
    \item{model}{Indicates the model the test uses: 
Use \code{model = "logistic"} for a two-valued outcome \code{Y} 
(the default) or \code{model = "linear"} for a continuous 
outcome. If \code{model} is not supplied, globaltest will try to
determine the model from \code{Y}.}
    \item{levels}{If \code{Y} is a factor (or a category in the PhenoData slot of \code{X}) 
and contains more than 2 levels: \code{levels} is a vector of levels of \code{Y} to test. If 
\code{levels} is length 2: test these 2 groups against each other. 
If levels is length 1: test that level against the others.}
    \item{adjust}{Confounders or risk factors for which the test must 
be adjusted. Must be either a data frame or the names of 
covariates in the phenoData slot of the exprSet \code{X}} 
    \item{permutation}{A logical flag. If \code{TRUE} \code{nperm} 
permutations are used to calculate the p-value instead of the 
asymptotic formulas. Recommended for small sample size. Not 
possible if an adjusted globaltest is used.}
    \item{nperm}{The number permutations used. The default is 10,000. 
If a number is specified for \code{nperm}, \code{permutation} is 
automatically set to \code{TRUE}.}
    \item{sampling}{A logical flag. If \code{TRUE} \code{ndraws} 
random sets of genes are drawn with the same number of genes as 
the tested group. Using this draws, an extra column of output 
\code{comparative.p} is generated, reporting how many of these 
random sets have a lower p-value than the tested group.}
    \item{ndraws}{The number of random groups of genes to be drawn.
The default is 1,000. If a number is specified for \code{ndraws},
\code{sampling} is automatically set to \code{TRUE}.}
    \item{verbose}{Prints some progress information if set to \code{TRUE}}.}

\details{The Global Test tests whether a group of genes (of any 
size from one single gene to all genes on the array) is 
significantly associated with a clinical outcome. The group could 
be for example a known pathway, an area on the genome or the set 
of all genes. The test investigates whether samples with similar 
clinical outcomes tend to have similar gene expression patterns. 
For a significant result it is not necessary that the genes in 
the group have similar expression patterns, only that many of 
them are correlated with the outcome.}

\note{If the number of rows of a matrix \code{X} does not match 
the length of the vector \code{Y}, but the number of columns 
does, the matrix \code{X} given is tacitly replaced by 
\code{t(X)} to make \code{X} and \code{Y} match. A warning is 
printed if \code{X} is square.}

\value{The function returns an object of class 
\code{\link{gt.result-class}}.}

\references{J. J. Goeman, S. A. van de Geer, F. de Kort and J. C. 
van Houwelingen, 2004, \emph{A global test for groups of genes: 
testing association with a clinical outcome}, 
\emph{Bioinformatics} 20 (1) 93--99. See also the vignette 
Globaltest.pdf included with this package.}

\author{Jelle Goeman: \email{j.j.goeman@lumc.nl}; Jan Oosting}

\seealso{\code{\link{geneplot}}, \code{\link{sampleplot}}, 
\code{\link{permutations}}, \code{\link{checkerboard}}, 
\code{\link{regressionplot}}.}

\examples{
    data(exampleX)      # Expression data (40 samples; 1000 genes)
    data(exampleY)      # Clinical outcome for the 40 samples
    pathway1 <- 1:25    # A pathway contains genes 1 to 25
    pathway2 <- 26:50   # another pathway
    gt <- globaltest(exampleX, exampleY, list(pathway1,pathway2))
    gt
}
\keyword{htest}
