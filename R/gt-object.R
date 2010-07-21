#==========================================================
# CLASS DEFINITION *** CLASS DEFINITION *** CLASS DEFINITION
#==========================================================
setClassUnion("matrixOrNULL", c("matrix", "NULL"))
setClassUnion("listOrNULL", c("list", "NULL"))
setClassUnion("data.frameOrNULL", c("data.frame", "NULL"))


setClass("gt.object", 
  representation(
    result = "matrix",
    extra = "data.frameOrNULL",
    call = "call", 
    functions = "environment",#"list",
    subsets = "listOrNULL",
    structure = "listOrNULL",
    weights = "listOrNULL",
    alternative = "matrixOrNULL",
    null = "matrixOrNULL",
    directional = "logical",
    legend = "list",
    model = "character"
  ),
  prototype = list(
    extra = NULL,
    subsets = NULL,
    structure = NULL,
    weights = NULL,
    alternative = NULL,
    null = NULL,
    directional = FALSE
  )
)

#==========================================================
# Function "show" prints a "gt.object" object
#==========================================================
setMethod("show", "gt.object", function(object)
{
  if (!is.null(alias(object))) {
    chars <- nchar(alias(object))
    mp <- gt.options()$max.print
    alias(object) <- ifelse(chars > mp,
      paste(substr(alias(object),1,mp-3), "...", sep=""),
      alias(object))
  }
  print(result(object), digits = 3)
})

setGeneric("summary")
setMethod("summary", "gt.object", function(object, ...)
{
  cat("\"gt.object\" object from package globaltest\n\n")
  cat("Call:\n")
  cat(deparse(object@call), "\n\n")
  df <- object@functions$df()
  cat("Model:", object@model, "regression.\n")
  cat("Degrees of freedom:", df[1], "total;", df[2], "null;", df[2], "+", df[3], "alternative.\n")
  if ((!is.null(object@subsets)) || !is.null(object@weights))
    cat("Tested:", length(object), "subsets.\n") 
  nperms = object@functions$nperms()
  cat("Null distibution: ")
  if (nperms[1]) {
    cat(if(!nperms[2]) "all", nperms[1], if(nperms[2]) "random", "permutations.\n")
  } else {
    cat("asymptotic.\n")
  }
  cat("\n")
  show(object)
})


#==========================================================
# Functions to extract relevant information from 
# a gt.object object
#==========================================================
setGeneric("result", function(object, ...) standardGeneric("result"))
setMethod("result", "gt.object",
  function(object) { 
    if (is.null(object@extra))
      data.frame(object@result, check.names=FALSE)
    else
      data.frame(object@extra, object@result, check.names=FALSE)
})

setGeneric(".result", function(object, ...) standardGeneric(".result"))
setMethod(".result", "gt.object",
  function(object) 
    object@result
)

setGeneric("weights") 
setMethod("weights", "gt.object", function(object) {
  
  # base weights based on covariate variance
  X <- object@functions$getX()
  weights <- colSums(X*X)
  if (object@model == "multinomial")
    weights <- rowSums(matrix(weights, object@functions$df()[3]))
  names(weights) <- object@functions$cov.names()
                                  
  # find weights for specific weights and subsets chosen
  if (length(object@subsets) > 0) {
    weights <- lapply(object@subsets, function(set) weights[set])
  }
  if (length(object@weights) > 0) 
    if (is.list(weights)) {
      weights <- lapply(as.list(1:length(weights)), function(i)
        weights[[i]] * object@weights[[i]])
      names(weights) <- names(object@weights)  
    } else
      weights <- lapply(object@weights, function(wts) weights * wts)
                  
  # set the max weight to 1  
  if (is.list(weights))
    weights <- lapply(weights, function(wts) wts / max(wts))
  else
    weights <- weights / max(weights)
      
  # reduce a list of length 1 to a vector    
  if (is.list(weights) && length(weights) == 1)
    weights <- weights[[1]]

  weights      
})

#==========================================================
setGeneric("subsets", function(object, ...) standardGeneric("subsets"))
setMethod("subsets", "gt.object", function(object, ...) {
  object@subsets
})


#==========================================================
setGeneric("p.value", function(object, ...) standardGeneric("p.value"))
setMethod("p.value", "gt.object",
  function(object) {
    out <- .result(object)[,1,drop=FALSE]
    colnames(out) <- NULL
    drop(out)
  }
)

#==========================================================
setGeneric("z.score", function(object, ...) standardGeneric("z.score"))
setMethod("z.score", "gt.object",
  function(object) {
    out <- .result(object)
    colnames(out) <- NULL            
    (out[,2] - out[,3]) / out[,4]
  }
)


#==========================================================
setGeneric("size", function(object, ...) standardGeneric("size"))
setMethod("size", "gt.object",
  function(object) {
    .result(object)[,"#Cov"]
  }
)


#==========================================================
# The subsetting methods for "gt.object"
#==========================================================
setMethod("[", "gt.object", 
            function(x, i, j,...,drop) 
{
  if (all(i %in% names(x)) || 
          all(i %in% 1:length(x)) ||
          all(i %in% -1:-length(x)) ||
          (is.logical(i) && (length(i) == length(x)))) {
    x@extra <- x@extra[i, ,drop=FALSE]
    if (!is.null(x@subsets)) x@subsets <- x@subsets[i]
    if (!is.null(x@weights)) x@weights <- x@weights[i]
    if (!is.null(x@structure)) {
      selector <- 1:length(x)
      names(selector) <- names(x)
      sub <- names(selector[i])
      if (length(sub) < length(x)) {
        x@structure$ancestors <- x@structure$ancestors[sub]
        x@structure$ancestors <- lapply(x@structure$ancestors, function(xx) intersect(xx, sub))
        x@structure$ancestors <- x@structure$ancestors[sapply(x@structure$ancestors, length) > 0]
        x@structure$offspring <- x@structure$offspring[sub]
        x@structure$offspring <- lapply(x@structure$offspring, function(xx) intersect(xx, sub))      
        x@structure$offspring <- x@structure$offspring[sapply(x@structure$offspring, length) > 0]
      }
    }
    x@result <- x@result[i, ,drop = FALSE] 
    x
  } else {
    stop("invalid index set", call. = FALSE)
  }
})            

setMethod("[[", "gt.object", 
            function(x, i, j,...,exact) 
{
   x[i]
})

#==========================================================
# The length method for "gt.object"
#==========================================================
setMethod("length", "gt.object", 
            function(x) 
{
  nrow(.result(x))
})            

#==========================================================
# The length method for "gt.object"
#==========================================================
leafNodes <- function(object, alpha=0.05, type = c("focuslevel","inheritance")) {

  if (is.null(object@structure$ancestors)) {     # Infer from sets
    sets <- object@subsets
    ancestors <- new.env(hash=TRUE)
    offspring <- new.env(hash=TRUE)
    for (i in 1:length(sets)) {
      namei <- names(sets)[i]
      for (j in 1:length(sets)) {
        namej <- names(sets)[j]
        if (i != j && length(sets[[i]]) <= length(sets[[j]]) && all(sets[[i]] %in% sets[[j]])) {
          ancestors[[namei]] <- c(ancestors[[namei]], namej)
          offspring[[namej]] <- c(offspring[[namej]], namei)
        }
      }
    }
    object@structure$ancestors <- as.list(ancestors)
    object@structure$offspring <- as.list(offspring)
  }

  if (missing(type)) 
    type <- c("focuslevel","inheritance")
  else
    type <- match.arg(type)  
  type <- type[type %in% names(object@extra)]
  if (length(type) < 1)
    stop("no focus level or inheritance p-values in object.")
  if (length(type) > 1)
    stop("both focus level and inheritance p-values in object. Please specify type.")

  ps <- object@extra[[type]]
  significant <- names(object)[ps <= alpha]
  sign.anc <- unique(unlist(object@structure$ancestors[significant]))
  leaves <- setdiff(significant, sign.anc)

  object[leaves]
}


#==========================================================
# The names and alias methods for "gt.object" 
# (applies to pathwaynames)
#==========================================================
setMethod("names", "gt.object", 
            function(x) 
{
  rownames(.result(x))
})      


setMethod("names<-", "gt.object", 
            function(x, value) 
{
  rownames(x@result) <- value
  if (!is.null(x@extra)) rownames(x@extra) <- value
  if (!is.null(x@weights)) names(x@weights) <- value
  if (!is.null(x@subsets)) names(x@subsets) <- value
  x
})            

setGeneric("alias")
setMethod("alias", "gt.object",
  function(object) {
    object@extra$alias  
})

setGeneric("alias<-", function(x, value) standardGeneric("alias<-"))
setMethod("alias<-", "gt.object",
  function(x, value) {
    if (length(value) != length(x))
      stop("alias has length ", length(value), " but object has length ", length(x))
    extra <- x@extra
    extra$alias <- value
    x@extra <- as.data.frame(extra, stringsAsFactors=FALSE)
    if (is.null(rownames(x@extra))) rownames(x@extra) <- names(x)
    x
})



#==========================================================
# A sort method for "gt.object"
#==========================================================
setMethod("sort", matchSignature(signature(x = "gt.object"), sort),
  function(x, decreasing = FALSE ) {
    mch <- match(c("focuslevel", "inheritance"), names(x@extra))
    if (any(!is.na(mch))) {
      mch <- mch[!is.na(mch)][1]
      adjp <- x@extra[[mch]]
      ix <- order(adjp, p.value(x), -z.score(x), decreasing=decreasing)
    } else
      ix <- order(p.value(x), -z.score(x), decreasing=decreasing)
    x[ix]
  }
)

#==========================================================
# Model.matrix extracts the model matrix (only if x=TRUE)
#==========================================================
setMethod("model.matrix", matchSignature(signature(object = "gt.object"), model.matrix),
  function(object, ... ) {
    list(alternative = object@alternative, null = object@null)
  }
)


#==========================================================
# Multiple testing correction for "gt.object" object
#==========================================================
setGeneric("p.adjust")
setMethod("p.adjust", matchSignature(signature(p = "gt.object"), p.adjust),
  function(p, method = p.adjust.methods, n = length(p)) {
    method <- method[1]
    method <- p.adjust.methods[grep(method, p.adjust.methods, ign=T)]
    if(length(method)==(0))   # this is just to get a good error message
      method <- match.arg(method)
    if (is.null(p@extra))
      p@extra <- data.frame(matrix(,length(p),0), row.names=names(p))
    if (missing(n))
      p@extra[[method]] <- p.adjust(p.value(p), method=method)
    else
      p@extra[[method]] <- p.adjust(p.value(p), method=method, n=n)
    p
  }
)



#==========================================================
# Histogram method to visualize permutations
#==========================================================
setGeneric("hist", function(x,...) standardGeneric("hist"))
setMethod("hist", "gt.object", function(x, ...)  {

  gt.hist <- function(x, breaks=20, main="", xlab = "Permutation test statistics", ...) {

    if (length(x) > 1)
    stop("length(object) > 1. Please reduce to a single test result")
  
    if (is.null(x@weights)) 
      weights <- rep(1, size(x))
    else
      weights <- x@weights[[1]]
    if (is.null(x@subsets))
      subset <- seq_len(size(x))
    else
      subset <- x@subsets[[1]]
  
    recalculate <- x@functions$permutations(subset, weights)
    Qs <- recalculate$permS
    Q <- recalculate$S
    nperm <- length(Qs)
    hst <- hist(Qs, xlim = c(1.1 * min(0, Qs, Q), 1.1 * max(Qs, Q)), breaks = breaks, 
      main = main, xlab = xlab, ...)
    h <- max(hst$counts)
    arrows( Q, h/2, Q, 0 , lwd=2)
    text( Q, h/2, 'Observed\ntest\nstatistic' , pos=3)
  
    # No output
    invisible(list(statistic = Q, histogram = hst))
  }
  
  gt.hist(x,...)
})        


#==========================================================
# Graph plot for focus level and inheritance procedures
#==========================================================
draw <- function(object, alpha=0.05, type = c("focuslevel","inheritance"), names=FALSE, sign.only = FALSE, interactive = FALSE) {

  # check availablity of packages
  require("Rgraphviz") || stop("package \"Rgraphviz\" is not available.")

  # find ancestors and offspring if missing
  if (is.null(object@structure$ancestors)) {     # Infer from sets
    sets <- object@subsets
    ancestors <- new.env(hash=TRUE)
    offspring <- new.env(hash=TRUE)
    for (i in 1:length(sets)) {
      namei <- names(sets)[i]
      for (j in 1:length(sets)) {
        namej <- names(sets)[j]
        if (i != j && length(sets[[i]]) <= length(sets[[j]]) && all(sets[[i]] %in% sets[[j]])) {
          ancestors[[namei]] <- c(ancestors[[namei]], namej)
          offspring[[namej]] <- c(offspring[[namej]], namei)
        }
      }
    }
    object@structure$ancestors <- as.list(ancestors)
    object@structure$offspring <- as.list(offspring)
  }

  # find type if missing
  if (missing(type)) 
    type <- c("focuslevel","inheritance")
  else
    type <- match.arg(type)  
  type <- type[type %in% names(object@extra)]
  if (length(type) < 1)
    stop("no focus level or inheritance p-values in object.")
  if (length(type) > 1)
    stop("both focus level and inheritance p-values in object. Please specify type.")

  ps <- object@extra[[type]]
  significant <- names(object)[ps <= alpha]

  parents <- ancestors2parents(object@structure$ancestors)

  # make the graph object
  graph <- as.matrix(sapply(names(object), function(node) names(object) %in% parents[[node]]))
  rownames(graph) <- colnames(graph) <- names(object)
  if (sign.only) {
    graph <- graph[significant, significant,drop=FALSE]
    if (length(significant)==0)
      stop("no significant nodes to plot.")
  }
  graph <- as(graph, "graphNEL")
  
  nodes <- buildNodeList(graph)
  edges <- buildEdgeList(graph)
  nAttrs <- list()
  eAttrs <- list()

  # color significant nodes
  if (!sign.only) {
    signode <- sapply(nodes, name) %in% significant
    names(signode) <- names(nodes)
    sigedge <- (sapply(edges, from) %in% significant) & (sapply(edges, to) %in% significant)
    names(sigedge) <- names(edges)
    nodecolor <- ifelse(signode, "black", "#BBBBBB")
    nAttrs$color <- nodecolor
    nAttrs$fontcolor <- nodecolor
    
    edgecolor <- ifelse(sigedge, "black", "#BBBBBB")
    eAttrs$color <- edgecolor
  }
  
  # if no names, give the plot numbers in order
  # this requires plotting the graph twice
  if (!names) {
    nAttrs$label <- 1:length(names(nodes))
    names(nAttrs$label) <- names(nodes)
    pg <- agopen(graph, name="pg", nodeAttrs = nAttrs, edgeAttrs = eAttrs)
    x <- getNodeXY(pg)$x
    y <- getNodeXY(pg)$y
    ordering <- sort.list(order(-y, x))
    nAttrs$label <- ordering
    names(nAttrs$label) <- names(nodes)
    plot(graph, attrs = list(node=list(shape="rectangle")), nodeAttrs = nAttrs, edgeAttrs = eAttrs)
  } else
    plot(graph, attrs = list(node=list(shape="rectangle")), nodeAttrs = nAttrs, edgeAttrs = eAttrs)
    
  # Make the plot interactive if asked
  if (interactive) {
    cat("Click in the plot to see name and alias. Press escape to return.\n")
    flush.console()
    repeat {
      p <- locator(n = 1) 
      if (is.null(p))
        break()
      pg <- plot(graph, attrs = list(node=list(shape="rectangle")), nodeAttrs = nAttrs, edgeAttrs = eAttrs)
      x <- getNodeXY(pg)$x
      y <- getNodeXY(pg)$y
      distance <- abs(p$x - x) + abs(p$y - y)
      idx <- which.min(distance)
      legend("topleft", legend = paste(nAttrs$label[idx], names(object)[idx], alias(object)[idx]), bg = "white", box.lty=0)
    }
  }
  
  # return a legend if the plot has numbers
  if (!names) {
    legend <- cbind(name=names(object), alias=alias(object))[sort.list(ordering),]
    legend <- data.frame(legend, stringsAsFactors=FALSE)
  } else
    legend <- NULL
    
  invisible(legend)
}
