gtGO <- function(response, exprs, ..., id, annotation, probe2entrez, ontology = c("BP", "CC", "MF"),
                  minsize=1, maxsize=Inf,
                  multtest = c("Holm", "focuslevel", "BH", "BY"), focuslevel = 10,
                  sort = TRUE) {

  # some input checking on response (unavoidable?)
  call <- match.call()
  if (is(exprs, "ExpressionSet")) {
    data <- pData(exprs)
  } else {
    data <- eval(call$data)
  }
  if (is.matrix(data))  
    data <- data.frame(data)
  formula.response <- (length(call$response) > 1) && (call$response[[1]] == "~")
  if (!formula.response) {
    name.response <- deparse(call$response)
    response <- eval(call$response, data, parent.frame())
  }  

  # get the right annotation package
  if (missing(annotation)) {
    if (is(exprs, "eSet")) {
      annotation <- NULL
      annotation <- annotation(exprs)
    } else
      stop("argument \"annotation\" is missing with no default.")
  }
  if (substr(annotation, nchar(annotation)-2, nchar(annotation)) == ".db")
    annotation <- substr(annotation, 1, nchar(annotation)-3)
  package <- paste(annotation, ".db", sep="")
  require(package, character.only=TRUE) || stop("package ", package, " is not available")

  # check whether "org" package is given
  if (substr(annotation,1,4) == "org.")
    extension <- "GO2ALLEGS"
  else
    extension <- "GO2ALLPROBES"
  GOOBJECT <- eval(as.name(paste(annotation, extension, sep="")))

  # reduce the terms
  if (missing(id))
    id <- mappedkeys(GOOBJECT)
  myGOTERM <- lookUp(id, "GO", "TERM", load=TRUE)

  # get the right ontology/ies
  if (!missing(ontology)) {
    choose <- sapply(myGOTERM, Ontology) %in% ontology
    id <- id[choose]
    myGOTERM <- myGOTERM[choose]
  }

  # retrieve sets
  sets <- lookUp(id, annotation, extension)
  sets <- lapply(sets, function(st) if (all(is.na(st))) character(0) else st)
                                                         
  # map back
  if (!missing(probe2entrez)) {
    if (is.environment(probe2entrez) || is(probe2entrez, "AnnDbBimap")) probe2entrez <- as.list(probe2entrez)
    if (is.list(probe2entrez)) probe2entrez <- unlist(probe2entrez)
    sets <- lapply(sets, function(set) {
      names(probe2entrez)[probe2entrez %in% set]
    })
  }

  # reduce to covariates available in exprs
  if (is(exprs, "eSet"))
    probes <- featureNames(exprs)
  else
    if (gt.options()$transpose)
      probes <- rownames(exprs)
    else
      probes <- colnames(exprs)
  sets <- lapply(sets, function(set) intersect(set, probes))

  # size restrictions
  size <- sapply(sets, length)
  choose <- size >= minsize & size <= maxsize
  sets <- sets[choose]
  myGOTERM <- myGOTERM[choose]
  if (length(sets) == 0) stop("No GO terms wih size between \"minsize\" and \"maxsize\"")
                   
  # perform tests and do multiple testing
  if (length(sets) > 1) {
    multtest <- match.arg(multtest)
    if (multtest == "focuslevel") { 
      ancestors <- lapply(as.list(ontology), function(ont) {
        ext <- paste(ont, "ANCESTOR", sep="")
        GOOBJ <- eval(as.name(paste("GO", ont, "ANCESTOR", sep="")))
        ontid <- intersect(keys(GOOBJ), id)
        if (length(ontid) > 0) lookUp(ontid, "GO", ext) else list()
      })
      ancestors <- do.call(c, ancestors)
      offspring <- lapply(as.list(ontology), function(ont) {
        ext <- paste(ont, "OFFSPRING", sep="")
        GOOBJ <- eval(as.name(paste("GO", ont, "OFFSPRING", sep="")))
        ontid <- intersect(keys(GOOBJ), id)
        if (length(ontid) > 0) lookUp(ontid, "GO", ext) else list()
      })
      offspring <- do.call(c, offspring)      
      offspring <- sapply(offspring, function(os) if (all(is.na(os))) character(0) else os)
      if (is.numeric(focuslevel))
        focuslevel <- findFocus(sets, ancestors = ancestors, offspring = offspring, maxsize = focuslevel)
      res <- gt(response, exprs, ...)
      res <- focusLevel(res, sets=sets, focus=focuslevel, ancestors = ancestors, offspring = offspring)
    } else {
      res <- p.adjust(gt(response, exprs, ..., subsets = sets), method = multtest)
    }
  } else 
    res <- gt(response, exprs, ..., subsets = sets) 
                                          
  # add names
  alias(res) <- sapply(myGOTERM, function(mgt) if (is(mgt, "GOTerms")) Term(mgt) else "")
  alias(res)[is.na(alias(res))] <- ""

  if (sort)
    res <- sort(res)
  
  res
}


gtKEGG <- function(response, exprs, ..., id, annotation, probe2entrez, 
                  multtest = c("Holm", "BH", "BY"), sort = TRUE) {

  # some input checking on response (unavoidable?)
  call <- match.call()
  if (is(exprs, "ExpressionSet")) {
    data <- pData(exprs)
  } else {
    data <- eval(call$data)
  }
  if (is.matrix(data))  
    data <- data.frame(data)
  formula.response <- (length(call$response) > 1) && (call$response[[1]] == "~")
  if (!formula.response) {
    name.response <- deparse(call$response)
    response <- eval(call$response, data, parent.frame())
  }  

  # get the right annotation package
  if (missing(annotation)) {
    if (is(exprs, "eSet")) {
      annotation <- NULL
      annotation <- annotation(exprs)
    } else
      stop("argument \"annotation\" is missing with no default.")
  }
  if (substr(annotation, nchar(annotation)-2, nchar(annotation)) == ".db")
    annotation <- substr(annotation, 1, nchar(annotation)-3)
  package <- paste(annotation, ".db", sep="")
  require(package, character.only=TRUE) || stop("package ", package, " is not available")

  # check whether "org" package is given
  if (substr(annotation,1,4) == "org.")
    extension <- "PATH2EG"
  else
    extension <- "PATH2PROBE"
  KEGGOBJECT <- eval(as.name(paste(annotation, extension, sep="")))
                                   
  # default terms
  if (missing(id))   
    id <- mappedkeys(KEGGOBJECT)

  # retrieve sets
  sets <- lookUp(id, annotation, extension)

  # map back
  if (!missing(probe2entrez)) {
    if (is.environment(probe2entrez) || is(probe2entrez, "AnnDbBimap")) probe2entrez <- as.list(probe2entrez)
    if (is.list(probe2entrez)) probe2entrez <- unlist(probe2entrez)
    sets <- lapply(sets, function(set) {
      names(probe2entrez)[probe2entrez %in% set]
    })
  }

  # reduce to covariates available in exprs
  if (is(exprs, "eSet"))
    probes <- featureNames(exprs)
  else
    if (gt.options()$transpose)
      probes <- rownames(exprs)
    else
      probes <- colnames(exprs)
  sets <- lapply(sets, function(set) intersect(set, probes))

  # perform tests and do multiple testing
  if (length(sets) > 1) {
    multtest <- match.arg(multtest)
    res <- p.adjust(gt(response, exprs, ..., subsets = sets), method = multtest)
  } else 
    res <- gt(response, exprs, ..., subsets = sets) 
  
  # add names
  alias(res) <- unlist(lookUp(names(res), "KEGG", "PATHID2NAME", load=TRUE))
  alias(res)[is.na(alias(res))] <- ""

  if (sort)
    res <- sort(res)
  
  res
}

gtBroad <- function(response, exprs, ..., id, annotation, probe2entrez, collection,
                  category = c("c1", "c2", "c3", "c4", "c5"), 
                  multtest = c("Holm", "BH", "BY"), sort = TRUE) {

  # some input checking on response (unavoidable?)
  call <- match.call()
  if (is(exprs, "ExpressionSet")) {
    data <- pData(exprs)
  } else {
    data <- eval(call$data)
  }
  if (is.matrix(data))  
    data <- data.frame(data)
  formula.response <- (length(call$response) > 1) && (call$response[[1]] == "~")
  if (!formula.response) {
    name.response <- deparse(call$response)
    response <- eval(call$response, data, parent.frame())
  }  

  # get the right annotation package
  if (missing(annotation)) {
    if (is(exprs, "eSet")) {
      annotation <- NULL
      annotation <- annotation(exprs)
    } else
      stop("argument \"annotation\" is missing with no default.")
  }
  if (substr(annotation, nchar(annotation)-2, nchar(annotation)) == ".db")
    annotation <- substr(annotation, 1, nchar(annotation)-3)
  package <- paste(annotation, ".db", sep="")
  require(package, character.only=TRUE) || stop("package ", package, " is not available")

  # read the file
  if (missing(collection) || !is(collection, "GeneSetCollection"))
    stop("Please specify a \"collection\", created with the getBroadSets function")
  
  # Get the right categories
  if (!missing(category)) {
    pw.cat <- sapply(sapply(collection, GSEABase::collectionType), GSEABase::bcCategory)
    collection <- collection[pw.cat %in% category]
  }
  
  # Get the right sets
  if (!missing(id)) {
    collection <- collection[id]
  }

  # Map symbol identifiers to anotation-specific identifiers
  collection <- GSEABase::mapIdentifiers(collection, GSEABase::AnnotationIdentifier(annotation))
  sets <- lapply(collection, GSEABase::geneIds)
  names(sets) <- names(collection)

  # map to probe identifiers
  if (!missing(probe2entrez)) {
    if (is.environment(probe2entrez) || is(probe2entrez, "AnnDbBimap")) probe2entrez <- as.list(probe2entrez)
    if (is.list(probe2entrez)) probe2entrez <- unlist(probe2entrez)
    sets <- lapply(sets, function(st) {
      names(probe2entrez)[probe2entrez %in% st]
    })
  }

  # reduce to covariates available in exprs
  if (is(exprs, "eSet"))
    probes <- featureNames(exprs)
  else
    if (gt.options()$transpose)
      probes <- rownames(exprs)
    else
      probes <- colnames(exprs)
  sets <- lapply(sets, function(set) intersect(set, probes))

  # perform tests and do multiple testing
  if (length(sets) > 1) {
    multtest <- match.arg(multtest)
    res <- p.adjust(gt(response, exprs, ..., subsets = sets), method = multtest)
  } else 
    res <- gt(response, exprs, ..., subsets = sets) 
  
  if (sort)
    res <- sort(res)
  
  res
}


###################################################
# Function for using Anni concept profiles
###################################################
gtConcept <- function(response, exprs, ..., annotation, probe2entrez,
                  conceptmatrix, concept2name = "conceptID2name.txt",
                  entrez2concept = "entrezGeneToConceptID.txt", threshold = 1e-4,
                  share = TRUE, multtest = c("Holm", "BH", "BY"),
                  sort = TRUE) {

  # some input checking on response (unavoidable?)
  call <- match.call()
  if (is(exprs, "eSet")) {
    data <- pData(exprs)
  } else {
    data <- eval(call$data)
  }
  if (is.matrix(data))
    data <- data.frame(data)
  formula.response <- (length(call$response) > 1) && (call$response[[1]] == "~")
  if (!formula.response) {
    name.response <- deparse(call$response)
    response <- eval(call$response, data, parent.frame())
  }

  # get the right annotation package
  if (missing(annotation)) {
    if (is(exprs, "eSet")) {
      annotation <- NULL
      annotation <- annotation(exprs)
    } else
      stop("argument \"annotation\" is missing with no default.")
  }
  if (substr(annotation, nchar(annotation)-2, nchar(annotation)) == ".db")
    annotation <- substr(annotation, 1, nchar(annotation)-3)
  package <- paste(annotation, ".db", sep="")
  require(package, character.only=TRUE) || stop("package ", package, " is not available")

  # read the big concept profiles matrix
  nconcepts <- length(scan(conceptmatrix, what="", sep="\t", nlines=1, quiet=TRUE))
  conceptmatrix <- scan(conceptmatrix, what="", sep="\t", quiet=TRUE)
  conceptmatrix <- matrix(conceptmatrix, length(conceptmatrix)/nconcepts, nconcepts, byrow=TRUE)
  rnames <- conceptmatrix[-1,1]
  cnames <- conceptmatrix[1,-1]
  conceptmatrix <- apply(conceptmatrix[-1,-1],2,as.numeric)
  rownames(conceptmatrix) <- rnames
  colnames(conceptmatrix) <- cnames

  # use the threshold
  conceptmatrix[conceptmatrix < threshold] <- 0

  # prepare entrez2concept
  entrez2concept <- matrix(scan(entrez2concept, what="", sep="\t", skip=1, quiet=TRUE), ncol=2, byrow=TRUE)
  nms <- entrez2concept[,1]
  entrez2concept <- as.list(entrez2concept[,2])
  names(entrez2concept) <- nms
                    
  # prepare concept2probe mapping
  if (is(exprs, "eSet"))
    probes <- featureNames(exprs)
  else
    probes <- rownames(exprs)
  if (missing(probe2entrez)) {
    if (substr(annotation,1,3) == "org") {
      probe2entrez <- as.list(probes)    # assuming probe ids are entrez ids
      names(probe2entrez) <- probes
    } else
      probe2entrez <- lookUp(probes, annotation, "ENTREZID")
  } else
    if (is.environment(probe2entrez) || is(probe2entrez, "AnnDbBimap"))
      probe2entrez <- as.list(probe2entrez)
  probe2entrez <- lapply(probe2entrez, as.character)
  probe2concept <- entrez2concept[unlist(probe2entrez)]
  names(probe2concept) <- names(probe2entrez)
  probe2concept <- probe2concept[sapply(probe2concept, length) > 0]

  # find duplicity of concepts and correct weights
  concept2probe <- turnListAround(probe2concept)
  duplicity <- sapply(concept2probe[rnames], length)
  names(duplicity) <- rnames
  if (share)     # divides weight of concept equally over multiple probes
    conceptmatrix <- conceptmatrix / matrix(duplicity, nrow = nrow(conceptmatrix), ncol=ncol(conceptmatrix))
  conceptmatrix <- conceptmatrix[duplicity>0,]

  # find weights
  weights <- lapply(probe2concept, function(pb)
    if(pb %in% rownames(conceptmatrix))
      conceptmatrix[pb,]
    else
      NULL)
  weights <- do.call(rbind, weights)
  weights <- apply(weights, 2, function(cp) cp[cp>0])

  # find alias
  concept2name <- matrix(scan(concept2name, what="", sep="\t", skip=1, quote="", quiet=TRUE), ncol=2, byrow=TRUE)
  nms <- concept2name[,1]
  concept2name <- as.list(concept2name[,2])
  names(concept2name) <- nms
                           browser()
  # perform tests and do multiple testing
  if (length(weights) > 1) {
    multtest <- match.arg(multtest)
    res <- p.adjust(gt(response, exprs, ..., weights = weights), method = multtest)
  } else
    res <- gt(response, exprs, ..., weights = weights)

  # add names
  alias(res) <- concept2name[names(res)]

  if (sort)
    res <- sort(res)

  res
}


