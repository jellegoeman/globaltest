setClass("GOstructure", representation(
  ids = "character",
  ancestors  = "list",
  offspring = "list",
  genesets = "list")
)

setMethod("show", "GOstructure", function(object) {
  cat("A GO structure object for globaltest with", length(object), "GO terms.\n")
})

setMethod("genesets", "GOstructure", function(object, ...) {
  object@genesets
})

setMethod("length", "GOstructure", function(x) {
  length(x@ids)
})



makeGOstructure <- function(data, annotation, top, only.ids, ontology = c("BP", "CC", "MF"), unreliable) {

  ll2name<-function(ids){
    ids<-ids[ids %in% annotation]
    nms<-data[match(ids,annotation)]
    names(nms)<-names(ids)
    nms
  }
  GO2ALLPROBES <-function(goid) {
    # return named vector of rownames
    if (length(annotation)==1) {
      lookUp(goid, annotation, "GO2ALLPROBES")
    } else {
      ll <- mget(goid, GOALLENTREZID)
      lapply(ll,ll2name)
    }
  }
  # Get the right ontology
  ontology <-  match.arg(ontology)
  if (missing(top)) {
    top <- switch(ontology, 
      "BP" = "GO:0008150",
      "CC" = "GO:0005575", 
      "MF" = "GO:0003674")
  }
  if (class(data) %in% c("exprSet","ExpressionSet")) {
    data<-featureNames(data)
  } else if (!is.character(data)) {
    stop("Invalid input of data")
  }
  if (length(annotation)>1 && length(annotation)!=length(data)) stop("Invalid annotation")

  ONancestor  <- paste(ontology, "ANCESTOR", sep="")
  ONoffspring  <- paste(ontology, "OFFSPRING", sep="")  
  require(GO)
  annotaton<-as.character(annotation)
  # Get the genes and remove unreliable annotations
  allgenes <- GO2ALLPROBES(top)[[1]]
  if (!missing(unreliable)) {
    allgenes <- allgenes[!is.na(names(allgenes))]
    unreliable <- match.arg(unreliable, c("IC","IDA","IEA","IEP","IGI","IMP","IPI","ISS","NAS","ND","RCA","TAS","NR"), several.ok = TRUE)
    allgenes <- allgenes[!(names(allgenes) %in% unreliable)]
  }
  allgenes <- intersect(allgenes, data)

  # Get all GO IDs
  ids <- lookUp(top, "GO", ONoffspring)
  if (is.list(ids)) ids<-ids[[1]]
  ids <- c(top, ids)
  if (!missing(only.ids)) {
    ids <- intersect(ids, only.ids)
  }
  # Find the sets and remove terms with empty sets

  genesets <- GO2ALLPROBES(ids)
  genesets <- lapply(genesets, function(set) intersect(set, allgenes))
  setsize <- sapply(genesets, length)
  ids <- ids[setsize > 0]
  genesets <- genesets[setsize > 0]
  if (length(ids) == 0)
    stop("Empty GOstructure")
  # Remove Synonyms
  idTERMs <- lookUp(ids, "GO", "TERM")
  allsynonyms <- unlist(sapply(idTERMs, Secondary))
  ids <- ids[!(ids %in% allsynonyms)]

  # Find ancestors and offspring
  ancestors <- lookUp(ids, "GO", ONancestor)
  ancestors <- lapply(ancestors, function(anc) setdiff(anc, "all"))
  offspring <- lookUp(ids, "GO", ONoffspring)
  offspring <- lapply(offspring, function(off) intersect(off, ids))
  
  # Make the output
  gg <- new("GOstructure")
  gg@ids <- ids
  gg@genesets <- genesets
  gg@ancestors <- ancestors
  gg@offspring <- offspring
  gg
}


getFocus <- function(GOstructure, maxatoms = 10) {

  # Find all nodes with at most 5*maxatoms offspring sets
  # and at most maxatoms atoms
  allnodes <- GOstructure@ids
  offspring <- GOstructure@offspring
  doable <- sapply(offspring, length) <= 5 * maxatoms
  for (term in allnodes[doable]) {
    if (length(offspring[[term]]) > 0) {
      atoms <- getAtoms(GOstructure@genesets[offspring[[term]]])
    } else {
      atoms <- character(0)
    }
    doable[term] <- (length(atoms) <= maxatoms)
  }
  
  # Find all doable nodes with no doable ancestors
  alldoable <- allnodes[doable]
  ancestors <- GOstructure@ancestors[doable]
  redundant <- sapply(alldoable, function(id) {
    any(ancestors[[id]] %in% alldoable)
  })
  alldoable[!redundant]
}


# This function extracts the smallest possible number of atoms from 
# sets, such that all sets can be built as unions of these atoms
# Note: These atoms lead to slightly conservative closed testing
# example: 1,2,23 is reduced to 1,2,3 and leads to unnessesary testing of 13 prior to testing 1.
# But: 1,12,13,23 is now reduced to 1,2,3 while otherwise it would not be reduced.
getAtoms <- function(allsets) {
  if (!is(allsets, "list")) {                                           # Coerce to a list
    allsets <- list(allsets)
  }
  allsets <- lapply(allsets, function(set) set[!is.na(set)])            # Remove NAs
  allsets <- allsets[sapply(allsets, length) > 0]                       # Remove empty sets
  setdiffs <- vector("list", length(allsets))                           # Prepare output
  i <- 0
  lastAtom <- character(0)
  while (length(allsets) > 0) {
    allsets <- lapply(allsets, function(set) setdiff(set, lastAtom))    # Remove probes in the last atom from every set
    allsets <- allsets[sapply(allsets, length) > 0]                     # Remove empty sets
    if (length(allsets) > 0) {
      i <- i+1
      allsets <- allsets[sort.list(sapply(allsets, length))]            # Sort by size
      lastAtom <- allsets[[1]]                                          # Next atom: smallest non-empty set
      setdiffs[i] <- allsets[1]                                         # Store the atom
    }
  } 
  setdiffs[1:i]
}
