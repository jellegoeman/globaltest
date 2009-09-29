\name{Multiple testing on the GO graph}

  \alias{makeGOstructure}
  \alias{getFocus}


\title{Multiple testing on the GO graph (Deprecated)}

\description{Three function to test (part of) the GO graph for
association of the gene expression profile of GO terms with a
response variable. Used together, these functions return
multiplicity-adjusted p-values calculated using the Focus Level
procedure that preserves the structure of the GO graph.}


\usage{
  makeGOstructure(data, annotation, top, only.ids, ontology = c("BP", "CC", "MF"), entrez, unreliable)
  getFocus(GOstructure, maxatoms = 10)}

\arguments{
    \item{data}{The data set for which the \code{\link{GOstructure}} object is to be made. Can be either an
    \code{\link[Biobase:class.ExpressionSet]{ExpressionSet}} or a vector of
    probe indicators.}
    \item{annotation}{The name of the metaData annotation package to be used. If a metadata package is not available, the name of the organism package (e.g. \code{org.HS.eg.db} can be given instead, and the mapping from probe identifiers to entrez can be supplied with the \code{entrez} argument.}
    \item{top}{The node to be used as top node of the \code{\link{GOstructure}}. Defaults
    to one of "GO:0008150", "GO:0005575", "GO:0003674", depending on he ontology chosen. The
    \code{\link{GOstructure}} object will only contain offspring of the chosen top node.}
    \item{only.ids}{A vector of GO ids. If this is supplied, GO ids not appearing in this list
    will not appear in the \code{\link{GOstructure}} object.}
    \item{ontology}{The ontology to be used. One of Biological Process (BP), Cellular
    Component (CC) or Molecular Function (MF).}
    \item{entrez}{A vector of EntrezGene ids with the same length as the number of genes in \code{data}, (and in the same order) giving the mapping from probe identifiers to EntrezGene ids. If this argument is used, the \code{annotation} argument should give the name of the organism annotation package.}
    \item{unreliable}{Can be used to designate one or more of the GO evidence codes as unreliable.
    Value must be a vector containing one or more from
    "IC", "IDA", "IEA", "IEP", "IGI", "IMP", "IPI", "ISS", "NAS", "ND", "RCA", "TAS", "NR".}
    \item{GOstructure}{An object of class \code{\link{GOstructure}}.}
    \item{maxatoms}{A tuning parameter that governs the choice of the focus level. \code{getFocus} will
    choose those GO terms in the focus level in such a way that all offspring gene sets of each focus level
    term can be constructed as unions of no more than \code{maxatoms} atom gene sets. Default value of
    \code{maxatoms} is 10. Higher values quickly lead to slower
    computation. Lower values typically lead to reduced power.}
}

\details{These functions should be used in the following order.
First use \code{makeGOstructure} to make a GO graph tailored to a
specific data set. Then \code{getFocus} can be used to choose a
focus level. Finally \code{gtGO} performs the focus level
procedure.}

\note{\code{gtGO} cannot be used in combination with the permutation
version of \code{\link{globaltest}}.}

\value{The function returns a named vector of multiplicity-adjusted
p-values. Adjusted p-values of GO terms not appearing in this vector
are larger than the chosen value of \code{maxalpha}.}

\references{For references, type: \code{citation("globaltest")}. See
also the vignette GlobalTest.pdf included with this package.}

\author{Jelle Goeman: \email{j.j.goeman@lumc.nl}; Jan Oosting}

\seealso{Examples in the vignette! \code{\link{globaltest}},
\code{\link{GOstructure}}, \code{\link{gt.multtest}}.}


\keyword{htest}
