\name{GOstructure}

  \alias{GOstructure}
  \alias{GOstructure-class}
  \alias{length,GOstructure-method}
  \alias{genesets,GOstructure-method}
  \alias{show,GOstructure-method}


\title{Class "GOstructure"  (Deprecated)}

\description{The class \code{GOstructure} stores the Gene Ontology
information that can be used by the \code{\link{gtGO}} function to
find a GO subgraph significantly associated with a response
variable. Objects of type \code{GOstructure} are typically made with
a call to \code{\link{makeGOstructure}}.}

\section{Slots}{
  \describe{
    \item{\code{ids}:}{Object of class "character". Vector of GO identifiers.}
    \item{\code{offspring}:}{A named list of vectors of GO identifiers. Lists all offspring terms of each term.}
    \item{\code{ancestors}:}{A named list of vectors of GO identifiers. Lists all ancestor terms of each term.}
    \item{\code{genesets}:}{A named list of vectors of gene identifiers. Lists the genes annotated to each GO id.}
    }
}

\section{Methods}{
  \describe{
    \item{show}{(GOstructure): Summarizes the object.}
    \item{length}{(GOstructure): Gives the number of GO ids in the object.}
    \item{genesets}{(GOstructure): Extracts the genesets slot.}
    }
}

\author{Jelle Goeman: \email{j.j.goeman@lumc.nl}; Jan Oosting}

\seealso{\code{\link{gtGO}}, \code{\link{makeGOstructure}}.}

\keyword{methods}
