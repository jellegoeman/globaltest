\name{gt gene set testing methods}

\alias{gtGO}
\alias{gtKEGG}
\alias{gtBroad}
\alias{gtConcept}

\title{Gene set testing of gene set databases using Global Test}

\description{A collection of procedures for performing the Global Test on gene set databases. Three function are provided for KEGG, for Gene Ontology and for the Broad Institute's gene sets.}

\usage{
gtKEGG (response, exprs, ..., id, annotation, probe2entrez,
           multtest = c("Holm", "BH", "BY"), sort = TRUE)

gtGO (response, exprs, ..., id, annotation, probe2entrez,
           ontology = c("BP", "CC", "MF"), minsize=1, maxsize=Inf,
           multtest = c("Holm", "focuslevel", "BH", "BY"),
           focuslevel = 10, sort = TRUE)


gtConcept (response, exprs, ..., annotation, probe2entrez,
           conceptmatrix, concept2name = "conceptID2name.txt",
           entrez2concept = "entrezGeneToConceptID.txt", threshold = 1e-4,
           share = TRUE, multtest = c("Holm", "BH", "BY"),
           sort = TRUE)

gtBroad (response, exprs, ..., id, annotation, probe2entrez, collection,
           category = c("c1", "c2", "c3", "c4", "c5"),
           multtest = c("Holm", "BH", "BY"), sort = TRUE)

}

\arguments{
    \item{response}{The response variable of the regression model. This is passed on to the \code{response} argument of \code{\link{gt}}.}
    \item{exprs}{The expression measurements. May be
      \code{\link[Biobase:class.ExpressionSet]{ExpressionSet}} or
      matrix. Passed on to the \code{alternative} argument of
      \code{\link{gt}}.}

    \item{...}{Any other arguments are also passed on to \code{\link{gt}}.}
    \item{id}{The identifier(s) of gene sets to be tested (character vector). If omitted, tests all gene sets in the database.}
    \item{annotation}{The name of the probe annotation package for the microarray that was used, or the name of the genome wide annotation package for the species (e.g. org.Hs.eg.db for human). If an organism package is given, the argument \code{probe2entrez} must be supplied. If \code{annotation} is missing, the function will  attempt to retrieve the annotation information from the \code{exprs} argument.}
    \item{probe2entrez}{Use only if no probe annotation package is available. A mapping from probe identifiers to entrez gene ids. May be an environment, named list or named vector.}
    \item{multtest}{The method of multiple testing correction. Choose from: Benjamini and Hochberg FDR control (BH); Benjamini and Yekutieli FDR control (BY) or Holm familywise error fontrol (Holm). For \code{gtGO} also the focus level method is available. See \code{\link{focusLevel}}.}
    \item{sort}{If \code{TRUE}, sorts the results to increasing p-values.}
    \item{ontology}{The ontology or ontologies to be used. Default is to use all three ontologies.}
    \item{minsize}{The minimum number of probes that may be annotated to a gene set. Gene sets with fewer annotated probes are discarded.}
    \item{maxsize}{The maximum number of probes that may be annotated to a gene set. Gene sets with more annotated probes are discarded.}
    \item{focuslevel}{The focus level to be used for the focus level method. Either a vector of gene set ids, or a numerical level. In the latter case, \code{\link{findFocus}} is called with \code{maxsize} at the specified level to find a focus level.}
    \item{collection}{The Broad gene set collection, created by a call to \code{\link[GSEABase:getObjects]{getBroadSets}}.}
    \item{conceptmatrix}{The name of the file containing the importance weights, i.e. concept profile associations between Anni concepts. In the matrix contained in the file, columns correspond to testable concepts, and rows correspond to entrez-concepts. Useable files can be downloaded from \url{http://www.biosemantics.org/weightedglobaltest}.}
    \item{concept2name}{The name of the file containing a mapping between Anni concepts and entrez identifiers. Useable files can be downloaded from \url{http://www.biosemantics.org/weightedglobaltest}.}
    \item{entrez2concept}{The name of the file containing a mapping between Anni concept numbers and names. Useable files can be downloaded from \url{http://www.biosemantics.org/weightedglobaltest}.}
    \item{threshold}{The relevance threshold for importance weights. Importance weights below the threshold are treated as zero.}
    \item{share}{If \code{TRUE}, the function divides the importance weight of a gene over all probes corresponding to the same entrez identifier. If \code{FALSE}, all probes get the full importance weight of the gene.}
    \item{category}{The subcategory of the Broad collection to be tested. The default is to test all sets.}
}

\details{These are utility functions to make it easier to do gene set testing of gene sets available in gene set databases. The functions automatically retrieve the gene sets, preprocess and select them, perform global test, do multiple testing correction, and sort the results on the basis of their p-values.

The four functions use different databases for testing. \code{gtKEGG} and \code{gtGO} use KEGG (\url{http://www.genome.jp/kegg}) and GO (\url{http://www.geneontology.org}); \code{gtConcept} uses the Anni database (\url{http://www.biosemantics.org/anni}), and \code{gtBroad} uses the MSigDB database (\url{http://http://www.broadinstitute.org/gsea/msigdb}). The \code{gtConcept} function differs from the other three in that it uses association weights between 0 and 1 for genes within sets, rather than having a hard cut-off for membership of a gene in a set.

All functions require that \code{annotate} and the appropriate annotation packages are installed. \code{gtKEGG} additionally requires the \code{KEGG.db} package; \code{gtGO} requires the \code{GO.db} package; \code{gtBroad} requires the user to download the XML file "msigdb_v2.5.xml"
from \\ \code{http://www.broad.mit.edu/gsea/downloads.jsp}, and to preprocess that file using the \code{\link[GSEABase:getObjects]{getBroadSets}} function. \code{gtConcept} requires files that can be downloaded from \url{http://www.biosemantics.org/weightedglobaltest}.
}

\value{A \code{\link{gt.object}}.}

\references{
Goeman, Van de Geer, De Kort and Van Houwelingen (2004). A global test for groups of genes: testing association with a clinical outcome. Bioinformatics 20 (1) 93-99.

Goeman, Oosting, Cleton-Jansen, Anninga and Van Houwelingen (2005). Testing association of a pathway with survival using gene expression data. Bioinformatics 21 (9) 1950-1957.
}

\author{Jelle Goeman: \email{j.j.goeman@lumc.nl}; Jan Oosting}

\seealso{
The \code{\link{gt}} function. The \code{\link{gt.object}} and useful functions associated with that object.
}

\examples{
    # Examples in the Vignette
}

\keyword{htest}
