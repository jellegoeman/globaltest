\name{globaltest}

\alias{globaltest}

\title{Global Test}

\description{In microarray data, tests whether a given group of 
genes is significantly associated with a given clinical variable.}

\usage{ globaltest(X, Y, test.genes = NULL, model = 'logistic', 
                            permutation = FALSE, nperm = NULL) }

\arguments{
    \item{X}{A matrix of gene expression data, where columns correspond to 
samples and rows to genes. The data should be properly normalized 
beforehand (and log- or otherwise transformed), but missing 
values are allowed (coded as \code{NA}). Gene and sample names 
can be included as the row and column names of \code{X}.}
    \item{Y}{A vector with the clinical outcome of interest, having one value
for each sample.} 
    \item{test.genes}{Indicates the group of genes to be tested. \code{test.genes} can 
be given in two formats. Either it can be a vector with 1 
(\code{TRUE}) or 0 (\code{FALSE}) for each gene in \code{X}, with 
1 indicating that the gene belongs to the group. Or it can be a 
vector containing the column numbers (in \code{X}) of the genes 
belonging to the group. The default is to test all genes.}
    \item{model}{Indicates the model the test uses: Use \code{model = 'logistic'} 
for a two-valued outcome \code{Y} (the default) or \code{model = 
'linear'} for a continuous outcome.}
    \item{permutation}{A logical flag. If \code{TRUE} permutations are used to 
calculate the p-value, if \code{FALSE} (the default) the 
formula's based on asymptotic theory are used. Choose 
\code{permutation = TRUE} for small sample sizes.}
    \item{nperm}{The number of permutations to be used. The default is 10,000. If 
a number is specified for \code{nperm}, \code{permutation} is 
automatically set to \code{TRUE}.} }

\details{The Global Test tests whether a group of genes (of any 
size from one single gene to all genes on the array) is 
significantly associated with a clinical outcome. The group could 
be for example a known pathway, an area on the genome or all 
genes. For a significant result it is not important if the genes 
in the group are positively or negatively correlated with the 
outcome. More precisely the test investigates whether samples 
with similar clinical outcomes tend to have similar gene 
expression patterns.}

\note{If the number of rows of \code{X} does not match the length 
of the vector \code{Y}, but the number of columns does, the 
matrix \code{X} given is tacitly replaced by \code{t(X)} to make 
\code{X} and \code{Y} match. A warning is printed if \code{X} is 
square.}

\value{The function returns an object of type \code{gt.result}.}

\references{J. J. Goeman, S.A. van de Geer, F. de Kort and J. C. 
van Houwelingen, \emph{A global test for association of a group 
of genes with a clinical outcome}, Technical Report MI 2003-03, 
Mathematical Institute, Leiden University. Available from 
\url{http://www.math.leidenuniv.nl/~jgoeman}.}

\author{Jelle Goeman: \email{j.j.goeman@lumc.nl}.}

\seealso{\code{\link{geneplot}}, \code{\link{regressionplot}}, 
\code{\link{checkerboard}}, \code{\link{permutations}}, 
\code{\link{gt.result}}.}

\examples{ 
    data(exampleX)      # Expression data (40 samples; 1000 genes)
    data(exampleY)      # Clinical outcome for the 40 samples
    pathway <- 1:25     # A pathway contains genes 1 to 25
    gt <- globaltest(exampleX, exampleY, test.genes = pathway)
    gt
    
    # For comparison: test all genes:
    gt <- globaltest(exampleX, exampleY)
    gt
    
    # For comparison: the permutation version:
    gt <- globaltest(exampleX, exampleY, test.genes = pathway, 
                                                permutation = TRUE)
    gt
    permutations(gt)
}

\keyword{htest}
    
  
    
