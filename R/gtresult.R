#==========================================================
# CLASS DEFINITION *** CLASS DEFINITION *** CLASS DEFINITION
#==========================================================
setClass("gt.result", 
  representation(
    res = "matrix", 
    eX = "matrix",
    genesets = "list",
    fit = "list", # A list of 1 item which is a "lm", c("lm", "glm") or "coxph" object
    method = "numeric",
    IminH = "matrix",
    PermQs = "matrix",
    samplingZs = "list"
  ),
  prototype = list(
    res = matrix(,0,0),
    eX = matrix(,0,0),
    pData = data.frame(),
    fit = list(NULL),
    method = 1,
    IminH = matrix(,0,0),
    PermQs = matrix(,0,0),
    samplingZs = list()
  )
)

#==========================================================
# PUBLIC METHODS *** PUBLIC METHODS *** PUBLIC METHODS
#==========================================================

#==========================================================
# Function "show" prints a "gt.result" object
#==========================================================
setMethod("show", "gt.result", function(object)
{
  npathways <- .nPathways(object)
  nsamples <- .nSamples(object)
  ngenes <- .nGenes(object)

  cat("Global Test result:\n")
  if (npathways == 1)
    cat("Data:", nsamples, "samples with", ngenes, "genes; 1 pathway tested\n")
  else
    cat("Data:", nsamples, "samples with", ngenes, "genes;", npathways, "pathways tested\n")
  if (!is.null(fit(object))) {
    cat("Model:", .model(object))
    if (.adjusted(object)) {
      cat(",", as.character(.formula(object))[c(2,1,3)], "\n")
      if (!is.na(.Rsquare(object)))
        cat("Adjusted:", 100 * round(.Rsquare(object), 3), "% of variance of Y remains after adjustment\n")
    } else
      cat("\n")
  }
  cat("Method:", .method(object), "\n")
  if (npathways > 0) {
    cat("\n") 
    print(signif(object@res, digits = 5))
  }
})


#==========================================================
# Functions to extract relevant information from 
# a gt.result object
#==========================================================
setMethod("result", "gt.result",
            function(object) object@res)


#==========================================================
setMethod("p.value", "gt.result",
            function(x) x@res[,6])


#==========================================================
setMethod("z.score", "gt.result",
            function(x) (x@res[,3] - x@res[,4]) / x@res[,5])


#==========================================================
# Function "fit" extracts the fitted adjustmodel
#==========================================================
setMethod("fit", "gt.result", 
            function(x) x@fit[[1]])



#==========================================================
# The subsetting method for "gt.result"
#==========================================================
setMethod("[", "gt.result", 
            function(x, i, j,...,drop) 
{
  if (all(i %in% names(x@genesets)) | 
          all(i %in% 1:length(x@genesets)) | 
          (is.logical(i) & (length(i) == length(x@genesets)))) {
    x@res <- x@res[i, ,drop = FALSE] 
    x@genesets <- x@genesets[i]
    x@PermQs <- x@PermQs[i, ,drop = FALSE]
    x
  } else {
    stop("invalid geneset", call. = FALSE)
  }
})            


#==========================================================
# The length method for "gt.result"
#==========================================================
setMethod("length", "gt.result", 
            function(x) 
{
  length(x@genesets)
})            


#==========================================================
# The names methods for "gt.result" 
# (applies to pathwaynames)
#==========================================================
setMethod("names", "gt.result", 
            function(x) 
{
  names(x@genesets)
})      


#==========================================================
setMethod("names<-", "gt.result", 
            function(x, value) 
{
  names(x@genesets) <- value
  rownames(x@res) <- value
  rownames(x@PermQs) <- value
  x
})            


#==========================================================
# A sort method for "gt.result"
#==========================================================
setMethod("sort", matchSignature(signature(x = "gt.result"), sort),
  function(x, partial = NULL, na.last = TRUE, decreasing = FALSE, 
      method = c("shell", "quick", "radix"), index.return = FALSE) {
    ix <- sort.list(p.value(x), partial, na.last, decreasing, method)
    x <- x[ix]
    if (index.return) 
      list(x = x, ix = ix)
    else
      x
  }
)


#==========================================================
# Combines two or more gt.result objects to one
#==========================================================
setMethod("combine", signature(x = "gt.result"),
  function(x, y) {
    if (!missing(y)) {
      if (!.formula(x) == .formula(y))
        stop("Different models.")
      if (!all(dim(x@eX) == dim(y@eX)) || !all(x@eX == y@eX))
        stop("Different gene expression matrices.")
      if (ncol(x@PermQs) != ncol(y@PermQs))
        stop("Different number of permutations used.")
      overlap <- intersect(rownames(x@res), rownames(y@res))
      if (length(overlap) > 0) {
        if (any(result(x)[overlap,1:6] != result(y)[overlap,1:6]))
          stop("Conflict: different genesets with the same name.")
        y@res <- y@res[!rownames(y@res) %in% overlap,]
      }
      if (.hasComparativeP(x) && !.hasComparativeP(y)) {
        y@res <- cbind(y@res, NA)
        colnames(y@res)[ncol(y@res)] <- "Comparative p"
      } else if (!.hasComparativeP(x) && .hasComparativeP(y)) {
        x@res <- cbind(x@res, NA)
        colnames(x@res)[ncol(x@res)] <- "Comparative p"
      }
      x@res <- rbind(x@res, y@res)
      x@genesets <- c(x@genesets, y@genesets)
      x@PermQs <- rbind(x@PermQs, y@PermQs)
      sampled <- unique(c(names(x@samplingZs), names(y@samplingZs)))
      x@samplingZs <- lapply(as.list(sampled), function(n)
        c(x@samplingZs[[n]], y@samplingZs[[n]])
      )
    }
    x
  }
)


#==========================================================
# PRIVATE METHODS *** PRIVATE METHODS *** PRIVATE METHODS
#==========================================================


#==========================================================
.Y <- function(gt) {
  model = .model(gt)
  if (model == "linear") {
    Y <- fit(gt)$y - fitted.values(fit(gt))
  } else if (model == "logistic") {
    Y <- fit(gt)$y - fitted.values(fit(gt))
    if (any(Y[fit(gt)$model[,1] == .levels(gt)[1]] < 0)) {
      Y <- -Y
    }
  } else if (model == "survival") {
    n <- .nSamples(gt)
    times <- as.vector(fit(gt)$y)[1:n]
    d <- as.vector(fit(gt)$y)[(n+1):(2*n)] 
    expci <- exp(fit(gt)$linear.predictors)
    dtimes <- unique(times[d == 1]) 
    nd <- length(dtimes)
    matrixO <- outer(times, dtimes, "==") * matrix(d, n, nd)
    atrisk <- outer(times, dtimes, ">=")
    hazinc <- as.vector(1 / (expci %*% atrisk))
    matrixP <- outer(expci, hazinc, "*") * atrisk
    matrixPO <- matrixP %*% t(matrixO)
    Y <- rowSums(matrixPO) - d
  } else if (model == "multinomial") {
    Y <- fit(gt)@y - fitted.values(fit(gt))
  }
  Y
}

#==========================================================
.levels <- function(gt) {
  if (.model(gt) == "logistic") {
    out <- levels(fit(gt)$model[,1])
    if (is.null(out)) 
      out <- c(1,0)
  } else {
    out <- colnames(.Y(gt))
    if (is.null(out)) 
      out <- 1:ncol(.Y(gt))
  }
  out
}


#==========================================================
.method <- function(gt) {
  switch (gt@method, 
    "No method",
    "Gamma approximation",
    "Asymptotic distribution",
    paste("All", ncol(gt@PermQs), "permutations"),
    paste(ncol(gt@PermQs), "random permutations")
  )
}

#==========================================================
.adjusted <- function(gt) {
  nrow(gt@IminH) > 0
}


#==========================================================
.IminH <- function(gt) {
  gt@IminH
}


#==========================================================
.genesets <- function(gt) {
  gt@genesets
}

#==========================================================
.Q <- function(gt) {
  result(gt)[,3]
}

#==========================================================
.nTested <- function(gt) {
  result(gt)[,2]
}

#==========================================================
.model <- function(gt) {
  if (is(fit(gt), "glm")) {
    model <- "logistic"
  } else if (is(fit(gt), "lm")) {
    model <- "linear"
  } else if (is(fit(gt), "coxph") || is(fit(gt), "coxph.null")) {
    model <- "survival"
  } else if (is(fit(gt), "mlogit")) {
    model <- "multinomial"
  }
  model
}

#==========================================================
.formula <- function(gt) {
  model <- .model(gt)
  if (model %in% c("logistic", "survival")) {
    ff <- fit(gt)$formula
  } else if (model == "multinomial") {
    ff <- fit(gt)@formula
  } else if (model == "linear") {
    ff <- formula(fit(gt)$terms)
  }
  ff
}

#==========================================================
.Rsquare <- function(gt) {
  model = .model(gt)
  if (model == "survival") {
    out = NA
  } else if (model == "logistic") {
    oldmu <- mean(fit(gt)$y)
    mu2 <- fit(gt)$weights
    n <- length(mu2)
    out <- sum(mu2) / (n * oldmu * (1 - oldmu))
  } else if (model == "linear") {
    Y <- residuals(fit(gt))
    nullY <- fit(gt)$y - mean(fit(gt)$y)
    out <- sum(Y * Y) / sum(nullY * nullY)
  } else if (model == "multinomial") {
    oldmu <- colMeans(fit(gt)@y)
    oldmu2 <- oldmu * (1-oldmu)
    mu <- fitted.values(fit(gt))
    mu2 <- mu * (1-mu)
    n <- nrow(mu)
    out <- sum(mu2) / (n * sum(oldmu2))
  }
  out
}

#==========================================================
.nSamples <- function(gt) {
  ncol(gt@eX)
}

#==========================================================
.nGenes <- function(gt) {
  nrow(gt@eX)
}

#==========================================================
.nPathways <- function(gt) {
  nrow(gt@res)
}

#==========================================================
.hasComparativeP <- function(gt) {
  "Comparative p" %in% colnames(result(gt))
}
