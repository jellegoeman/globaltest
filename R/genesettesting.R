gtGO <- function(response, exprs, ..., id, annotation, probe2entrez, ontology = c("BP", "CC", "MF"),
                  minsize=1, maxsize=Inf,
                  multtest = c("Holm", "focuslevel", "BH", "BY"), focuslevel = 10,
                  sort = TRUE) {

  # some input checking on response (unavoidable?)
  call <- match.call()
  if (is(exprs, "ExpressionSet")) {
    require("Biobase") || stop("ExpressionSet input but Biobase package not available")
    data <- pData(exprs)
  } else {
    data <- eval(call$data)
  }
  if (is.matrix(data))  
    data <- data.frame(data)
  formula.response <- (length(call$response) > 1) && (call$response[[1]] == "~")
  if (!formula.response) {
    if (length(call$response) > 1 && deparse(call$response[1]) == "Surv()")
      require("survival") || stop("Surv input but survival package not available")
    name.response <- deparse(call$response)
    response <- eval(call$response, data, parent.frame())
  }  

  # get the right annotation package
  if (missing(annotation) && is(exprs, "eSet")) {
    annotation <- NULL
    annotation <- annotation(exprs)
  }
  if (substr(annotation, nchar(annotation)-2, nchar(annotation)) == ".db")
    annotation <- substr(annotation, 1, nchar(annotation)-3)
  package <- paste(annotation, ".db", sep="")
  require(package, character.only=TRUE) || stop("package ", package, " is not available")
  require(GO.db) || stop("package GO.db is not available")
  require(annotate) || stop("package annotate is not available")

  # check whether "org" package is given
  if (substr(annotation,1,4) == "org.")
    extension <- "GO2ALLEGS"
  else
    extension <- "GO2ALLPROBES"
  GOOBJECT <- eval(as.name(paste(annotation, extension, sep="")))

  # reduce the terms
  if (missing(id))
    id <- mappedkeys(GOOBJECT)
  myGOTERM <- lookUp(id, "GO", "TERM")

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
    probes <- rownames(exprs)
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
      ancestors <- unlist(sapply(ontology, function(ont) {
        ext <- paste(ont, "ANCESTOR", sep="")
        GOOBJ <- eval(as.name(paste("GO", ont, "ANCESTOR", sep="")))
        ontid <- intersect(keys(GOOBJ), id)
        if (length(ontid) > 0) lookUp(ontid, "GO", ext) else list()
      }), recursive=FALSE)
      names(ancestors) <- substring(names(ancestors), 4)
      offspring <- unlist(sapply(ontology, function(ont) {
        ext <- paste(ont, "OFFSPRING", sep="")
        GOOBJ <- eval(as.name(paste("GO", ont, "OFFSPRING", sep="")))
        ontid <- intersect(keys(GOOBJ), id)
        if (length(ontid) > 0) lookUp(ontid, "GO", ext) else list()
      }), recursive=FALSE)      
      names(offspring) <- substring(names(offspring), 4)
      offspring <- sapply(offspring, function(os) if (all(is.na(os))) character(0) else os)
      if (is.numeric(focuslevel))
        focuslevel <- findFocus(sets, ancestors = ancestors, offspring = offspring, maxsize = focuslevel)
      res <- gt(response, exprs, ...)
      res <- focusLevel(res, sets=sets, focus=focuslevel, ancestors = ancestors, offspring = offspring)
    } else {
      res <- multtest(gt(response, exprs, ..., subsets = sets), method = multtest)
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
    require("Biobase") || stop("ExpressionSet input but Biobase package not available")
    data <- pData(exprs)
  } else {
    data <- eval(call$data)
  }
  if (is.matrix(data))  
    data <- data.frame(data)
  formula.response <- (length(call$response) > 1) && (call$response[[1]] == "~")
  if (!formula.response) {
    if (length(call$response) > 1 && deparse(call$response[1]) == "Surv()")
      require("survival") || stop("Surv input but survival package not available")
    name.response <- deparse(call$response)
    response <- eval(call$response, data, parent.frame())
  }  

  # get the right annotation package
  if (missing(annotation) && is(exprs, "eSet")) {
    annotation <- NULL
    annotation <- annotation(exprs)
  }
  if (substr(annotation, nchar(annotation)-2, nchar(annotation)) == ".db")
    annotation <- substr(annotation, 1, nchar(annotation)-3)
  package <- paste(annotation, ".db", sep="")
  require(package, character.only=TRUE) || stop("package ", package, " is not available")
  require(KEGG.db) || stop("package KEGG.db is not available")
  require(annotate) || stop("package annotate is not available")

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
    probes <- rownames(exprs)
  sets <- lapply(sets, function(set) intersect(set, probes))

  # perform tests and do multiple testing
  if (length(sets) > 1) {
    multtest <- match.arg(multtest)
    res <- multtest(gt(response, exprs, ..., subsets = sets), method = multtest)
  } else 
    res <- gt(response, exprs, ..., subsets = sets) 
  
  # add names
  alias(res) <- unlist(lookUp(names(res), "KEGG", "PATHID2NAME"))
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
    require("Biobase") || stop("ExpressionSet input but Biobase package not available")
    data <- pData(exprs)
  } else {
    data <- eval(call$data)
  }
  if (is.matrix(data))  
    data <- data.frame(data)
  formula.response <- (length(call$response) > 1) && (call$response[[1]] == "~")
  if (!formula.response) {
    if (length(call$response) > 1 && deparse(call$response[1]) == "Surv()")
      require("survival") || stop("Surv input but survival package not available")
    name.response <- deparse(call$response)
    response <- eval(call$response, data, parent.frame())
  }  

  # get the right annotation package
  if (missing(annotation) && is(exprs, "eSet")) {
    annotation <- NULL
    annotation <- annotation(exprs)
  }
  if (substr(annotation, nchar(annotation)-2, nchar(annotation)) == ".db")
    annotation <- substr(annotation, 1, nchar(annotation)-3)
  package <- paste(annotation, ".db", sep="")
  require(package, character.only=TRUE) || stop("package ", package, " is not available")
  require(GSEABase) || stop("package GSEABase is not available")

  # read the file
  if (missing(collection) || ! is(collection, "GeneSetCollection"))
    stop("Please specify a \"collection\", created with the getBroadSets function")
  
  # Get the right categories
  if (!missing(category)) {
    pw.cat <- sapply(sapply(collection, collectionType), bcCategory)
    collection <- collection[pw.cat %in% category]
  }
  
  # Get the right sets
  if (!missing(id)) {
    collection <- collection[id]
  }

  # Map symbol identifiers to anotation-specific identifiers
  collection <- mapIdentifiers(collection, AnnotationIdentifier(annotation))
  sets <- lapply(collection, geneIds)
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
    probes <- rownames(exprs)
  sets <- lapply(sets, function(set) intersect(set, probes))

  # perform tests and do multiple testing
  if (length(sets) > 1) {
    multtest <- match.arg(multtest)
    res <- multtest(gt(response, exprs, ..., subsets = sets), method = multtest)
  } else 
    res <- gt(response, exprs, ..., subsets = sets) 
  
  if (sort)
    res <- sort(res)
  
  res
}
