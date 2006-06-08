\name{Multiple testing on the GO graph}

  \alias{makeGOstructure}
  \alias{getFocus}
  \alias{gtGO}


\title{Multiple testing on the GO graph}

\description{Three function to test (part of) the GO graph for
association of the gene expression profile of GO terms with a
response variable. Used together, these functions return
multiplicity-adjusted p-values calculated using the Focus Level
procedure that preserves the structure of the GO graph.}


  \usage{
  gtGO(..., GO, focus, maxalpha = 0.05, stopafter = 100, verbose = FALSE)
  makeGOstructure(data, annotation, top, only.ids, ontology = c("BP", "CC", "MF"), unreliable)
  getFocus(GOstructure, maxatoms = 10)}

\arguments{
    \item{...}{Arguments describing the tests to be preformed are passed on to \code{\link{globaltest}}.}
    \item{GO}{An object of class \code{\link{GOstructure}} describing the structure of the GO
    graph. This object should be created using \code{makeGOstructure}}
    \item{focus}{A vector of GO ids to describe te focus level. Typically made using \code{getFocus}.}
    \item{maxalpha}{The maximum multiplicity-adjusted p-value. The algorithm will stop when this value is exceeded.}
    \item{stopafter}{The maximum number of significant GO terms to be found. The algorithm will stop when this value is exceeded.}
    \item{verbose}{If set to \code{TRUE}, prints much more extensive progress information.}
    \item{data}{The data set for which the \code{\link{GOstructure}} object is to be made. Can be either an
    \code{\link[Biobase:exprSet-class]{exprSet}} or a vector of
    probe indicators.}
    \item{annotation}{The name of the metaData annotation package to be used.}
    \item{top}{The node to be used as top node of the \code{\link{GOstructure}}. Defaults
    to one of "GO:0008150", "GO:0005575", "GO:0003674", depending on he ontology chosen. The
    \code{\link{GOstructure}} object will only contain offspring of the chosen top node.}
    \item{only.ids}{A vector of GO ids. If this is supplied, GO ids not appearing in this list
    will not appear in the \code{\link{GOstructure}} object.}
    \item{ontology}{The ontology to be used. One of Biological Process (BP), Cellular
    Component (CC) or Molecular Function (MF).}
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

\examples{
    # Breast cancer data (exprSet) from the Netherlands Cancer
    # Institute with small GO structure (all descendants of Cell Cycle)
    # See the vignette for details on making a GO structure object.
    data(vandeVijver)
    data(gostructure.vandeVijver)

    focus <- getFocus(gostructure.vandeVijver)
    adjustedP <- gtGO(vandeVijver, "Surv(TIMEsurvival, EVENTdeath)",
      GO = gostructure.vandeVijver, focus = focus)

    # The significant GO graph can be visualized using
    # the GOstats and Rgraphviz packages.
    # See the vignette for details. Type: vignette("GlobalTest")
}


\keyword{htest}
