#==========================================================
# Function "globaltest" performs the global test on (a list of) 
#    subsets of the data
#
# X is a data matrix (n rows are samples, p columns are genes) 
#    or Biobase exprSet.
# NB: if the dimension of X does not fit the dimension of Y, 
#    but t(X) does, t(X) is used.
# Code missing values as NA
#
# Y is a vector of n clinical outcomes or the name or index of one of 
#    the phenoData variables ( if X is an exprSet) 
#
# test.genes = list of vectors
#       vector can be a length p vector of 0 and 1
#           (1 = selected; 0 = not selected).
#       or a vector with the numbers of the selected genes
#           e.g. c(5,8,9) to test genes 5, 8 and 9.
#       or a vector with the ids in the rownames of X
#           e.g c("AA173143","AA461092","AA487442","AA490243","R26706","R61845")
# 
# OPTIONS:
# model = 'logistic': for two-value Y.
# model = 'linear': for continuous Y.
# model = 'survival'.
# d = event indicator for the survival model
# event = value of d which indicates the event
# levels: vector of groups to test in Y. Not needed if Y is binominal
#           If Y contains > 2 levels then the following methods are used
#           if levels contains 1 value: this groups is tested against all other samples
#           if levels contains 2 values: these groups are tested against each other
# adjust = data.frame: confounders or riskfactors for which the test must be adjusted
#       may be coded as names or indices of variables in the PhenoData slot of X (if exprSet)
#
# RESULT
# array with 7 columns containing
# p.value, Q, EQ, seQ, comparative.p, length of vector, rows found in X
# where comparative.p is the fraction of random pathways of the same size 
#   with a lower p-value
#==========================================================

globaltest <- function(X, Y, test.genes,
                        model, levels,
                        d, event = 1,
                        adjust, ...)


{
  # 0: the default for adjust:
  if (missing(adjust))
    adjust <- Y ~ 1
    
  # 1: make matrix X and pDataX:
  if (is(X, "exprSet")) {
    require("Biobase")
    if (is.data.frame(adjust)) {
      if (any(sampleNames(X) != rownames(adjust)))
        stop("samplenames of X and adjust inconsistent", call. = FALSE)
      if (any(names(adjust) %in% names(pData(X))))
        stop("duplicate variable names in pData(X) and adjust", call. = FALSE)
      pDataX <- cbind(adjust, pData(X))
      X <- exprs(X)
      if (is.null(colnames(X)))
        colnames(X) <- rownames(pDataX)
      adjust <- names(adjust)
    } else {
      pDataX <- pData(X)
      rownames(pDataX) <- sampleNames(X)
      X <- exprs(X)
      colnames(X) <- rownames(pDataX)
    }
  } else { # X is not an exprSet
    if (!is.vector(X) & !is.data.frame(X) & !is.matrix(X))
      stop("X should be of type 'matrix'.", call. = FALSE)
    X <- as.matrix(X)
    if (is.data.frame(adjust)) {
      pDataX <- adjust
      adjust <- names(adjust)
    } else {
      pDataX <- NULL
    }
  }
    
  # 2: coerce Y into a vector
  if ( is(Y, "formula") ) {
    adjust <- Y
    Y <- as.character(Y[2])
  }
  if ( length(Y)==1 ) {
    if (Y %in% names(pDataX)) {
      nameY <- Y
      Y <- pDataX[[Y]]
    } else { 
      stop("Y not among pData(X) variables", call. = FALSE)
    }
  } else {
    nameY <- NULL
    if ( !is.vector(Y) & !is.matrix(Y) & !is.data.frame(Y) & !is.factor(Y) )
        stop("Y should be of type 'vector'", call. = FALSE)
    if ( is.data.frame(Y) )
      Y <- as.matrix(Y)
    if ( is.matrix(Y) ) {
      if ( nrow(Y) == 1 )
        Y <- t(Y)
      if ( ncol(Y) == 1 ) {
        nameY <- colnames(Y)
        samplenamesY <- rownames(Y)
        Y <- as.vector(Y)
        names(Y) <- samplenamesY
      } else {
        stop("Y should be a single vector.", call. = FALSE)
      }
    }
    if (is.null(nameY))
      nameY <- "Y" 
  }
  
  # 3: Determine the model from the input Y
  if (missing(model)) {
    if (is.factor(Y) | (length(unique(Y)) <= 2) | !missing(levels)) {
      model <- 'logistic'
    } else {
      if (is.numeric(Y)) {
        if (!missing(d)) 
          model <- 'survival'
        else
          model <- 'linear'            
      } else {
        stop("model could not be determined from input Y: please specify model", call. = FALSE)
      }
    } 
  }
  model <- match.arg(model, c('linear','logistic','survival'))
  if (model == 'survival')
    require("survival")

  # 4: Preparation of X and Y in the logistic model using option levels
  if (model=='logistic') {
    if (missing(levels)){
      # Only 2 levels should be here, test for 1 now, later checks will find other errors
      levels <- NULL
      levels <- levels(factor(Y))
      if (length(levels)==1) 
        stop("There should be more than 1 group in the data.", call. = FALSE)
    } else {
      if (!all(levels %in% levels(factor(Y))))
        stop("input 'levels' incompatible with levels of Y", call. = FALSE)
      }
    if (length(levels)==2 ) {
    # create a subset of samples
      X <- as.matrix(X[,Y == levels[1] | Y == levels[2]])
      pDataX <- pDataX[Y == levels[1] | Y == levels[2], ,drop = FALSE]
      Y <- Y[Y == levels[1] | Y == levels[2]]
    }
    if (length(levels)<=2 ) {
      samplenamesY <- names(Y)
      Y <- (Y != levels[1])
      names(Y) <- samplenamesY
    }else{ 
      stop("No more than 2 groups can be tested. Use option: levels.", call. = FALSE)    
    }
    Y <- 0 + Y
  } else {
    levels <- numeric(0)
  }

  # 5: check if X and Y are numeric
  if (!is.numeric(X) | !is.numeric(Y))
      stop("X and Y must be numeric.", call. = FALSE)
  if (any(is.na(Y)))
      stop("missing values not allowed in Y.", call. = FALSE)
  
  # 6: check dimensions of X, Y and pDataX
  n <- length(Y)
  if (n == 1)
    stop("only one sample", call. = FALSE)
  if (ncol(X) == n) {
      # probably all exprsets will be transposed here
      p <- nrow(X)
      X <- t(X)
  }else{
      if (nrow(X) == n)
           p <- ncol(X)
      else
          stop("Dimensions of X and Y don't match.", call. = FALSE)
  }
  if (n == p & is.null(pDataX))
      warning("As many samples as genes.\n Columns of X are assumed samples; rows assumed genes.", call.=FALSE)
  if (is.null(pDataX)) {
    pDataX <- data.frame(Y)
    if ( (!is.null(rownames(X))) & (is.null( names(Y) )) ) {
      rownames(pDataX) <- rownames(X)
    }
  }
  if (nrow(pDataX) != n)
    stop("dimension of adjust incompatible with Y", call. = FALSE)

  # 7: extract genenames and compare samplenames
  genenames <- colnames(X)
  if (  (is.null(rownames(X))) & (!is.null(names(Y))) ) {
    rownames(X) <- names(Y)
  }else{
    if ( (!is.null(rownames(X))) & (is.null( names(Y) )) ) {
      names(Y) <- rownames(X)
    }else{ 
      if ( !all(names(Y) == rownames(X)) )
        stop("Sample names in X inconsistent with sample names in Y.", call. = FALSE)
    }
  }
  if (any(rownames(pDataX) != names(Y)))
    stop("samplenames in adjust inconsistent with samplenames in Y or X", call. = FALSE)

  # Check correct input in case of survival model
  if (model == 'survival') {
    named <- "d"
    if (missing(d)) {
      d <- rep(1, times = n)
    }else{
      if (!is.null(pDataX) & is.character(d) & (length(d) == 1) ) {
        if (d %in% names(pDataX)) {
          named <- as.character(d)
          d <- pDataX[,d]
        } else 
          stop("d not among pData(X) variables", call. = FALSE)
      }else{
        if (length(d) != n) {
          stop("length of d not equal to length of Y", call. = FALSE)
        }else{
          if (!all(names(d) == names(Y)))
            stop("names of d inconsistent with names of Y", call. = FALSE)  
          if (any(is.na(d)))
            stop("missing values not allowed in d", call. = FALSE)
        }
      }
      if (all(d %in% c(0,1)) & all(event == 1) ) {
        d <- as.numeric(d)
      } else {
        d <- as.numeric(d %in% event)
        named <- "d"
      }
      if (sum(d) == 0)
        stop("there are no events.", call. = FALSE)
    }
  }

  # Check for correct input of adjust
  if (is(adjust, "formula")) {
    if (!(as.character(adjust)[2] %in% c("Y", nameY)))
      stop("if a formula, adjust should begin with Y ~", call. = FALSE)
  } else {
    if (all(adjust %in% names(pDataX))) {
      adjust <- formula(paste("Y ~", paste(adjust, collapse = "+")))
    } else {
      if (any(adjust %in% names(pDataX)))
        stop("not all variables in adjust found among names(pData(X))", call. = FALSE)
      else
        stop("incorrect input of adjust", call. = FALSE)
    }
  }
  adjust <- formula(paste(nameY, "~", as.character(adjust)[3]))
  
  # remove redundant levels
  if (!is.null(pDataX)) {
    samples <- rownames(pDataX)
    pDataX <- data.frame(lapply(as.list(pDataX), function(x) { if (is.factor(x)) factor(x) else x}))
    rownames(pDataX) <- samples
  } 
  
  # Fit adjustmodel and adjust Y and X
  # Standardize Y in case of the linear model or unadjusted logistic (same variance for i=1,...,n)
  if (model == 'linear') {
    fit <- try(lm(adjust, pDataX, x = TRUE, na.action = 'na.fail'))
    if (is(fit, "try-error"))
        stop("Fitting the adjustmodel failed.", call.=FALSE)
    nullfit <- lm(Y ~ 1)
    Y <- residuals(fit)
    nullY <- residuals(nullfit)
    Z <- fit$x
    nadjust <- ncol(Z)
    if (df.residual(fit) > 0)
      mu2 <- sum(Y*Y)/df.residual(fit)
    else
      stop("no degrees of freedom left after adjustment", call. = FALSE)
    Rsquare <- sum(Y * Y) / sum(nullY * nullY)
    Y <- Y / sqrt(mu2)
    IminH <- diag(n) - Z %*% solve(t(Z) %*% Z, t(Z))
    parameters <- data.frame(matrix(,n,0))
  }
  if (model == 'logistic') {
    if (as.character(adjust[3]) != "1") {
      old.mu <- mean(Y)
      null.Y <- Y - old.mu
      fit <- try(glm(adjust, pDataX, family = binomial, x = TRUE, na.action = 'na.fail'))
      if (is(fit, "try-error"))
        stop("Fitting the adjustmodel failed.", call.=FALSE)   
      mu <- fitted.values(fit)
      Z <- fit$x
      nadjust <- ncol(Z)
      Y <- Y - mu
      mu2 <- mu * (1 - mu)
      Rsquare <- sum(mu2) / (n * old.mu * (1 - old.mu))
      IminH <- diag(n) - diag(mu2) %*% Z %*% solve(t(Z) %*% diag(mu2) %*% Z, t(Z))
      parameters <- data.frame(mu, mu2)
    } else {
      mu <- rep(mean(Y), times = n)
      mu2 <- mu * (1 - mu)
      Y <- (Y - mu) / sqrt(mu2)
      Rsquare <- 1
      Z <- matrix(rep(1, times = n))
      nadjust <- 1
      IminH <- diag(n) - matrix(rep(1/n, times = n*n), n, n)
      parameters <- data.frame(matrix(,n,0))
    }
  } 
  if (model == 'survival') {
    times <- Y
    adjust <- formula(paste("Surv(", nameY, ",", named, ") ~ ", as.character(adjust[3])))
    fit <- try(coxph(adjust, pDataX, x = TRUE, na.action = 'na.fail'))
    if (is(fit, "try-error"))
      stop("Fitting the adjustmodel failed.", call.=FALSE) 
    Z <- fit$x
    nadjust <- 1 + ncol(Z)
    expci <- exp(fit$linear.predictors)
    dtimes <- unique(times[d == 1]) 
    nd <- length(dtimes)
    matrixO <- outer(times, dtimes, "==") * matrix(d, n, nd)
    ties <- any(colSums(matrixO) > 1)
    atrisk <- outer(times, dtimes, ">=")
    hazinc <- as.vector(1 / (expci %*% atrisk))
    matrixP <- outer(expci, hazinc, "*") * atrisk
    matrixPO <- matrixP %*% t(matrixO)
    matrixM <- diag(d) %*% outer(times, dtimes, "<") - matrixPO %*% outer(times, dtimes, "<")
    matrixW <- diag(rowSums(matrixPO)) - matrixPO %*% t(matrixPO)
    parameters <- data.frame(times, d, expci)
    Y <- rowSums(matrixPO) - d
    if (ncol(Z) > 0)
      IminH <- diag(n) - matrixW %*% Z %*% solve(t(Z) %*% matrixW %*% Z, t(Z))
    else
      IminH <- diag(n)
    Rsquare <- as.numeric(NA)
  }
  #Rsquare <- as.numeric((null.Y %*% IminH %*% IminH %*% null.Y) / (null.Y %*% null.Y))
  if (!is.na(Rsquare) & (Rsquare < 0.01))
    stop("not enough variance in Y remaining after adjustment", call. = FALSE)
  
  # Impute missing values of X and rescale X
  col.meanX <- colMeans(X, na.rm = TRUE)
  X <- X - rep(1, times=n) %o% col.meanX
  X[is.na(X)] <- 0
  adjX <- t(IminH) %*% X
  if (model != "survival") {
    if ((model == "logistic") & (nadjust > 1)) 
      norm <- sqrt(sum(mu2 %*% (adjX * adjX)) / (10 * p))
    else {
      norm <- sqrt(sum(adjX * adjX) / (10 * p))
    }
    X <- X / norm
  }

  # Coerce test.genes into a list and check correct input
  if ( missing(test.genes) )
    test.genes <- list(1:p)
  
  if ( !is.list(test.genes) ) 
    test.genes <- list(test.genes)

  test.genes <- lapply(test.genes, function(tg) { 
    if ( !is.vector(tg) & !is.null(tg) )
      stop("test.genes should be a (list of) vector(s)", call. = FALSE)

    if ( all(tg %in% c(0,1)) & (length(tg) == p) ) {
      test.names <- names(tg)[tg == 1]
      if (  !is.null(genenames )  &  !is.null(test.names) ) {
        if( any( genenames[tg] != test.names ) ) 
          warning("Gene names in X inconsistent with gene names in test.genes.", call. = FALSE)
      }
      tg <- (1:p)[tg == 1]
    } else {
      if (is.character(tg))  {
        tg <- match(tg, genenames)
      } else {  
        if ( all(tg %in% 1:p) ) {
          test.names <- names(tg)
          names(tg) <- NULL
          if (  !is.null(genenames )  &  !is.null(test.names) ) {
            if( any( genenames[tg] != test.names ) ) 
            warning("Gene names in X inconsistent with gene names in test.genes.", call. = FALSE)
          }
        } else {
          stop("option 'test.genes' should be a (list of) vector(s) of gene names or numbers.", call. = FALSE)
        }
      } 
    }  
    tg
  })
  path.n <- sapply(test.genes, length) 
  test.genes <- lapply(test.genes, function(tg) if (!is.null(tg)) unique(tg[!is.na(tg)]) else tg ) 
  test.n <- sapply(test.genes, length)
  
  # Calculate the test statistic for each geneset
  res <- t(sapply(test.genes, function(tg) {
      
    out <- numeric(4)
    names(out) <- c("Q","EQ","seQ","p.val") 
    
    # select genes to be tested
    X.sel <- as.matrix(X[,tg])
    # Number of selected genes that are in array (test.genes can be larger than that (in res[index, "path.n"] earlier)
    m <- ncol(X.sel)
    if (m>0) {

      # calculate test statistic Q
      R <- ( t(IminH) %*% X.sel %*% t(X.sel) %*% (IminH) ) / m
      XX <- (X.sel %*% t(X.sel)) / m
      Q <- ( Y %*% XX %*% Y )
      if (model == "survival") {
        if (ties) {
          tiecorrect <- sum( sapply(1:length(dtimes), 
            function(i) {
              if (sum(matrixO[,i]) == 1)
                0
              else { 
                matrixtie <- outer(matrixO[,i], matrixO[,i]) - diag(matrixO[,i])
                sum(XX[as.logical(matrixtie)]) - 2 * sum(matrixP[,i] %*% XX %*% matrixtie) + sum(matrixtie) * (matrixP[,i] %*% XX %*% matrixP[,i])
              }
            }
          ))
          Q <- Q - tiecorrect
        }
      }
       
      # Expectation and variance of Q and p-value
      trR <- sum(diag(R))
      trRR <- sum(R*R)
      tr2R <- trR^2
      trR2 <- sum(diag(R*R))
      if (model == 'logistic') {
        if (as.character(adjust[3]) == "1") {
          EQ <- trR
          mu1 <- mu[1] 
          mumu <- mu2[1]
          K <- ( 1 - 6 * mu1 + 6 * mu1^2 ) / mumu
          varQ <- K * ( trR2 - tr2R / n ) + 2 * trRR - 2 * tr2R / (n-1)
        } else {
          # in this case Y is not standardized
          RV <- R %*% diag(mu2)
          EQ <- as.numeric(sum(diag(RV)))
          mu4 <- mu * (1-mu)^4 + (1-mu) * mu^4
          varQ <-  sum(diag(R) * diag(R) * (mu4 - 3 * mu2 * mu2)) + 2 * sum(diag(RV %*% RV))
        }
        seQ <- sqrt(varQ)
        scl <- varQ / (2 * EQ)
        dfr <- EQ / scl
        p.value <- pchisq(Q / scl, dfr, lower.tail = FALSE)
      }
      if (model == 'linear') {
        EQ <- sum(diag(R))
        varQ <- (2 / (n - nadjust + 2)) * ( (n - nadjust) * trRR - tr2R )
        seQ <- sqrt(varQ)
        scl <- varQ / (2 * EQ)
        dfr <- EQ / scl
        p.value <- pchisq(Q / scl, dfr, lower.tail = FALSE)
      }    
      if (model == 'survival') {
        EQ <- sum(diag(R %*% matrixW))
        tussen <- matrix(diag(R), n, ncol(matrixP)) + 2 * R %*% (matrixM - matrixP)
        matrixT <- tussen - matrix(colSums(matrixP * tussen), n, ncol(matrixP), byrow = TRUE)
        varQ <- sum(matrixPO * ((matrixT * matrixT) %*% t(matrixO)))
        Q <- (Q - EQ) / sqrt(varQ)
        p.value <- pnorm(-Q)
        EQ <- 0
        seQ <- 1
      }
      
      # the returns for this geneset
      out["p.val"] <- p.value
      out["Q"] <- Q
      out["EQ"] <- EQ
      out["seQ"] <- seQ
    }else{
      # returns for an empty geneset
      out["p.val"] <- NA
      out["Q"] <- NA
      out["EQ"] <- NA
      out["seQ"] <- NA
    }
  out  
  }))
  res <- cbind(path.n, test.n, res)
  rownames(res) <- names(test.genes)
  
  Qs <- matrix(,length(test.genes),0)
  rownames(Qs) <- names(test.genes)
  
  # returns an object of type "gt.result"
  gt <- new("gt.result",
    res = res, 
    X = X,
    Y = Y,
    test.genes = test.genes,
    IminH = IminH,
    pars = parameters,
    Rsquare = Rsquare,
    Qs = Qs,
    model = model,
    adjustmodel = adjust,
    levels = levels,
    df.adjust = as.integer(nadjust))
    
  # check deprecated options:
  dots <- list(...)
  sampling <- dots[["sampling"]]
  ndraws <- dots[["ndraws"]]
  permutation <- dots[["permutation"]]
  nperm <- dots[["nperm"]]
  if (!is.null(permutation) | !is.null(nperm)) {
    if (!is.null(nperm)) {
      gt <- permutation(gt, nperm = nperm)
    } else {
      gt <- permutation(gt)
    }
    warning("option permutation has been deprecated. Use function: permutations", call. = FALSE)
  }
  if (!is.null(sampling) | !is.null(ndraws)) {
    if (!is.null(ndraws)) {
      gt <- sampling(gt, ndraws = ndraws)
    } else {
      gt <- sampling(gt)
    }
    warning("option sampling has been deprecated. Use function: sampling", call. = FALSE)
  }
 
  gt
}
#==========================================================


#==========================================================
setClass("gt.result", representation(
    res = "matrix", 
    X = "matrix",
    Y = "vector",
    test.genes = "list",
    IminH = "matrix",
    pars = "data.frame",
    Rsquare = "numeric",
    Qs = "matrix",
    model = "character",
    adjustmodel = "formula",
    levels = "vector",
    df.adjust = "integer")
)
#==========================================================


#==========================================================
# Function "show" prints out a result of type "gt.result"
# such as results from a call to "globaltest"
#==========================================================

setMethod("show", "gt.result", function(object)
{
  npathways <- length(object@test.genes)
  nsamples <- length(object@Y)
  ngenes <- dim(object@X)[2]

  cat("Global Test result:\n")
  if (npathways == 1)
    cat("Data:", nsamples, "samples with", ngenes, "genes; 1 pathway tested\n")
  else
    cat("Data:", nsamples, "samples with", ngenes, "genes;", npathways, "pathways tested\n")
  cat("Model:", object@model)
  if (as.character(object@adjustmodel[3]) != "1") {
    cat(",", as.character(object@adjustmodel[2]), as.character(object@adjustmodel[1]), as.character(object@adjustmodel[3]), "\n")
    if (!is.na(object@Rsquare))
      cat("Adjusted:", 100 * round(object@Rsquare, 3), "% of variance of Y remains after adjustment\n")
  } else
    cat("\n")
  if (ncol(object@Qs) > 0)
    cat("Using", ncol(object@Qs), "permutations of Y\n")

  cat("\n")
    
  res <- data.frame(object@res)
  if ( ncol(res) < 7 ) {
    names(res) <- c("genes","tested","Statistic Q","Expected Q","sd of Q","p-value")
  }else{
    names(res) <- c("genes","tested","Statistic Q","Expected Q","sd of Q","p-value","comp. p")
  }
  print(signif(res, digits = 5))
})


#==========================================================
# Two functions to extract relevant information from 
# a gt.result object
#==========================================================
if( !isGeneric("result") )
    setGeneric("result", function(object) standardGeneric("result"))

setMethod("result", "gt.result",
            function(object) 
{
    object@res
})

#==========================================================
if( !isGeneric("p.value") )
    setGeneric("p.value", function(gt) standardGeneric("p.value"))

setMethod("p.value", "gt.result",
            function(gt) gt@res[,"p.val"])



#==========================================================
# The subsetting method for "gt.result"
#==========================================================
setMethod("[", "gt.result", 
            function(x, i, j,...,drop) 
{
  if (all(i %in% names(x@test.genes)) | 
          all(i %in% 1:length(x@test.genes)) | 
          (is.logical(i) & (length(i) == length(x@test.genes)))) {
    x@res <- x@res[i, ,drop = FALSE] 
    x@test.genes <- x@test.genes[i]
    x@Qs <- x@Qs[i, ,drop = FALSE]
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
  length(x@test.genes)
})            

#==========================================================
# The names methods for "gt.result" 
# (applies to pathwaynames)
#==========================================================
setMethod("names", "gt.result", 
            function(x) 
{
  names(x@test.genes)
})      

setMethod("names<-", "gt.result", 
            function(x, value) 
{
  names(x@test.genes) <- value
  rownames(x@res) <- value
  rownames(x@Qs) <- value
  x
})            

#==========================================================
# A sort method for "gt.result"
#==========================================================
if( !isGeneric("sort") ) setGeneric("sort")

setMethod("sort", matchSignature(signature(x = "gt.result", index.return = "logical"), sort),
  function(x, partial = NULL, na.last = NA, decreasing = FALSE, 
      method = c("shell", "quick", "radix"), index.return) {
    ix <- sort.list(p.value(x), partial, na.last, decreasing, method)
    x <- x[ix]
    if (index.return) 
      list(x = x, ix = ix)
    else
      x
  }
)
setMethod("sort", matchSignature(signature(x = "gt.result", index.return = "missing"), sort),
  function(x, partial = NULL, na.last = NA, decreasing = FALSE, 
      method = c("shell", "quick", "radix"), index.return) {
    sort(x, partial, na.last, decreasing, method, index.return = FALSE)
  }
)

  
#==========================================================
# Histogram method for permutations of "gt.result"
#==========================================================
if( !isGeneric("hist") )
    setGeneric("hist", function(x,...) standardGeneric("hist"))

setMethod("hist", "gt.result",
            function(x, ...) 
{
  gt.hist <- function(x, geneset = NULL, ...) {
    if (!is.null(geneset))
      x <- x[geneset]
    if (length(x@test.genes) > 1)
      stop("more than one geneset in x", call. = FALSE)
    if (ncol(x@Qs) == 0)
      stop("no permutations in x: try hist(permutations(x)) instead", call. = FALSE)
      
    Qs <- as.vector(x@Qs)
    Q <- x@res[,"Q"]
    nperm <- length(Qs)
    hst <- hist(Qs, xlim = c(1.1 * min(0, Qs, Q), 1.1 * max(Qs, Q)), breaks = 20, 
      main = "",
      xlab = "Values of Q for permuted Y", ...)
    h <- max(hst$counts)
    arrows( Q, h/5, Q, 0 )
    text( Q, h/4.5, 'Q' )

    # No output
    invisible(NULL)
  }
  gt.hist(x,...)
})            


#==========================================================
# Function "permutations" compares the theoretical values 
# of EQ, seQ and the p.value to values based on permutations 
# of the clinical outcome, which may be better for small 
# sample sizes. It summarizes the Q-values for the permutations 
# in a histogram, in which an arrow points out the value of the 
# true Q, for comparison
#==========================================================
permutations <- function(gt, geneset, nperm = 10^4)
{

    # check correct input of gt
    if ( !is(gt, "gt.result"))
      stop("permutations should be applied to a globaltest result", call. = FALSE)
    if (gt@df.adjust > 1)
      stop("the permutation procedure is not applicable for the adjusted global test", call. = FALSE)      
        
    # check correct input of nperm
    if ( !(nperm > 0) )
      stop("option nperm should be a positive integer", call. = FALSE)
    if (nperm <= ncol(gt@Qs))
      gt@Qs <-  gt@Qs[,1:nperm]
    else {
      nperm <- nperm - ncol(gt@Qs)
      # extract the right test.genes vector
      if (!missing(geneset))
        gt <- gt[geneset]
      QQs <- t(sapply(1:length(gt@test.genes), function(index) {
        if (gt@res[index, "test.n"] == 0)
          Qs <- rep(0, times = nperm)
        else {
          test.genes <- gt@test.genes[[index]]
          if (is.character(test.genes))
            test.genes <- intersect(test.genes, colnames(gt@X))
    
          # Recreate Y and XX
          X <- as.matrix(gt@X[,test.genes])
          m <- ncol(X)
          R <- ( t(gt@IminH) %*% X %*% t(X) %*% (gt@IminH) ) / m
          XX <- (X %*% t(X)) / m
          Y <- gt@Y
          n <- length(Y)
          Q <- gt@res[index, "Q"]
  
          if (gt@model == "survival") {
            # recreate matrices
            times <- gt@pars[["times"]]
            d <- gt@pars[["d"]]
            expci <- gt@pars[["expci"]]
            dtimes <- unique(times[d == 1]) 
            nd <- length(dtimes)
            matrixO <- outer(times, dtimes, "==") * matrix(d, n, nd)
            ties <- any(colSums(matrixO) > 1)
            atrisk <- outer(times, dtimes, ">=")
            hazinc <- as.vector(1 / (expci %*% atrisk))
            matrixP <- outer(expci, hazinc, "*") * atrisk
            matrixPO <- matrixP %*% t(matrixO)
            matrixM <- diag(d) %*% outer(times, dtimes, "<") - matrixPO %*% outer(times, dtimes, "<")
            matrixW <- diag(rowSums(matrixPO)) - matrixPO %*% t(matrixPO)
          }
  
          # Calculate Q for nperm permutations of Q
          if (gt@model != 'survival') {
            # Recalculate the Q-value for permutations of Y
            permQ <- function(npm) {
              pms <- apply( matrix(rnorm(n * npm), n, npm), 2, sort.list )
              Y.pm <-  matrix( Y[pms], n, npm )
              colSums(( XX %*% Y.pm ) * Y.pm)
            }

            chunk <- 5000
            nchunks <- trunc(nperm / chunk)
            rest <- nperm - nchunks * chunk
            if (rest == 0)
              chunks <- as.list(rep(chunk, nchunks))
            else
              chunks <- as.list(c(rep(chunk, nchunks), rest))
            Qs <- unlist(lapply(chunks, permQ))
          } else { 
            # survival model:
            if (!ties) {
              permQ <- function(npm) {
                pms <- apply( matrix(rnorm(n * npm), n, npm), 2, sort.list )
                Y.pm <-  matrix( Y[pms], n, npm )
                Q.pms <- colSums(( R %*% Y.pm ) * Y.pm)
                EQs <- rep(0, npm)
                varQs <- rep(0, npm)
                for (j in 1:nd) {
                  pj <- matrixP[,j]
                  PJ <- matrix( pj[pms], n, npm )
                  mj <- matrixM[,j]
                  MJ <- matrix( mj[pms], n, npm )
                  tussen <- matrix(diag(R), n, npm) + 2 * R %*% (MJ - PJ)
                  TJ <- tussen - matrix(colSums(PJ * tussen), n, npm, byrow = TRUE)
                  varQs <- varQs + colSums(PJ * TJ * TJ)
                  EQs <- EQs + diag(R) %*% PJ - colSums(( R %*% PJ ) * PJ)
                }
                (Q.pms - EQs) / sqrt(varQs)
              }
              chunk <- 5000
              nchunks <- trunc(nperm / chunk)
              rest <- nperm - nchunks * chunk
              if (rest == 0)
                chunks <- as.list(rep(chunk, nchunks))
              else
                chunks <- as.list(c(rep(chunk, nchunks), rest))
              Qs <- unlist(lapply(chunks, permQ))
            } else {
              # ties present: much slower generation of permutations
              Qs <- replicate(nperm, {
                shuffle <- sample(n)
                R <- R[shuffle, shuffle]
                XX <- XX[shuffle, shuffle]
                Q.sample <- as.numeric(Y %*% XX %*% Y)
                EQ.sample <- sum(R * matrixW)
                tiecorrect <- sum( sapply(1:length(dtimes), 
                  function(i) {
                    if (sum(matrixO[,i]) == 1)
                      0
                    else { 
                      matrixtie <- outer(matrixO[,i], matrixO[,i]) - diag(matrixO[,i])
                      sum(XX[as.logical(matrixtie)]) - 2 * sum(matrixP[,i] %*% XX %*% matrixtie) + sum(matrixtie) * (matrixP[,i] %*% XX %*% matrixP[,i])
                    }
                  }
                ))
                Q.sample <- Q.sample - tiecorrect
                tussen <- matrix(diag(R), n, ncol(matrixP)) + 2 * R %*% (matrixM - matrixP)
                matrixT <- tussen - matrix(colSums(matrixP * tussen), n, ncol(matrixP), byrow = TRUE)
                varQ <- sum(matrixPO * ((matrixT * matrixT) %*% t(matrixO)))
                (Q.sample-EQ.sample)/sqrt(varQ)
              })
            }
          }
        }
        Qs  
      })) 
      gt@Qs <- cbind(gt@Qs, QQs)
    }
    
    gt@res[,"EQ"] <- apply(gt@Qs, 1, mean)
    gt@res[,"seQ"] <- apply(gt@Qs, 1, sd)
    gt@res[,"p.val"] <- apply(gt@Qs >= matrix(gt@res[,"Q"], nrow(gt@Qs), ncol(gt@Qs)), 1, mean)
    gt@res <- gt@res[,1:6, drop = FALSE]
    
    gt
}
#==========================================================


#==========================================================
# Function "sampling" compares the p.value(s) found with 
# p-values from randomly generated "pathways" of the same size
# as the tested pathway.
#==========================================================
sampling <- function(gt, geneset, ndraws = 10^3)
{
  # check correct input of gt
  if ( !is(gt, "gt.result"))
    stop("sampling should be applied to a globaltest result", call. = FALSE)
  if (ncol(gt@Qs) > 0)
    stop("sampling cannot be applied to a permutation globaltest result", call. = FALSE)
      
  # check correct input of ndraws
  if ( !(ndraws > 0) )
    stop("option ndraws should be a positive integer", call. = FALSE)

  # extract the right test.genes vector
  if (!missing(geneset))
    gt <- gt[geneset]
   
  # recreate matrices
  p <- ncol(gt@X)
  n <- nrow(gt@X)
  if (gt@model == 'logistic') { 
     if (gt@df.adjust == 1) {  
       mu1 <- mean(gt@Y < 0)
       mumu <- mu1 * (1-mu1)
       K <- ( 1 - 6 * mu1 + 6 * mu1^2 ) / mumu 
     } else {
       mu <- gt@pars[["mu"]]
       mu2 <- gt@pars[["mu2"]]
       musq <- mu*mu
       mu4 <- mu - 4*musq + 6*musq*mu - 3*musq*musq  
     }
  }
  if (gt@model == "survival") {
    times <- gt@pars[["times"]]
    d <- gt@pars[["d"]]
    expci <- gt@pars[["expci"]]
    dtimes <- unique(times[d == 1]) 
    nd <- length(dtimes)
    matrixO <- outer(times, dtimes, "==") * matrix(d, n, nd)
    ties <- any(colSums(matrixO) > 1)
    atrisk <- outer(times, dtimes, ">=")
    hazinc <- as.vector(1 / (expci %*% atrisk))
    matrixP <- outer(expci, hazinc, "*") * atrisk
    matrixPO <- matrixP %*% t(matrixO)
    matrixM <- diag(d) %*% outer(times, dtimes, "<") - matrixPO %*% outer(times, dtimes, "<")
    matrixW <- diag(rowSums(matrixPO)) - matrixPO %*% t(matrixPO)
  }
     
  setsizes <- unique(sapply(gt@test.genes, length))
  setsizes <- setsizes[!(setsizes %in% c(0,p))]

  sampling.matrix <- t(sapply(setsizes, function(m) {
    replicate(ndraws, {
      if (m < p/2)
        X <- as.matrix(gt@X[,sample(1:p,m)])
      else
        X <- as.matrix(gt@X[,-sample(1:p,p-m)])
      R <- ( t(gt@IminH) %*% X %*% t(X) %*% (gt@IminH) ) / m
      XX <- (X %*% t(X)) / m
      Q <- as.numeric(gt@Y %*% XX %*% gt@Y)
      if (gt@model == 'logistic') {
        if (gt@df.adjust == 1) {
          trR <- sum(diag(R))
          trRR <- sum(R*R)
          tr2R <- trR*trR
          trR2 <- sum(diag(R*R))
          EQ <- trR
          varQ <- K * ( trR2 - tr2R / n ) + 2 * trRR - 2 * tr2R / (n-1)
        } else {
          RV <- R %*% diag(mu2)
          EQ <- sum(diag(RV))
          varQ <-  sum(diag(R) * diag(R) * (mu4 - 3 * mu2 * mu2)) + 2 * sum(diag(RV %*% RV))
        }
        scl <- varQ / (2 * EQ)
        dfr <- EQ / scl
        pout <- pchisq(Q / scl, dfr, lower.tail = FALSE)
      }
      if (gt@model == 'linear') {
        trR <- sum(diag(R))
        trRR <- sum(R*R)
        tr2R <- trR*trR
        EQ <- trR
        varQ <- (2 / (n - gt@df.adjust + 2)) * ( (n - gt@df.adjust) * trRR - tr2R )
        scl <- varQ / (2 * EQ)
        dfr <- EQ / scl
        pout <- pchisq(Q / scl, dfr, lower.tail = FALSE)
      }    
      if (gt@model == 'survival') {
        EQ <- sum(R * matrixW)
        tiecorrect <- sum( sapply(1:length(dtimes), 
          function(i) {
            if (sum(matrixO[,i]) == 1)
              0
            else { 
              matrixtie <- outer(matrixO[,i], matrixO[,i]) - diag(matrixO[,i])
              sum(XX[as.logical(matrixtie)]) - 2 * sum(matrixP[,i] %*% XX %*% matrixtie) + sum(matrixtie) * (matrixP[,i] %*% XX %*% matrixP[,i])
            }
          }
        ))
        Q <- Q - tiecorrect
        tussen <- matrix(diag(R), n, ncol(matrixP)) + 2 * R %*% (matrixM - matrixP)
        matrixT <- tussen - matrix(colSums(matrixP * tussen), n, ncol(matrixP), byrow = TRUE)
        varQ <- sum(matrixPO * ((matrixT * matrixT) %*% t(matrixO)))
        pout <- pnorm(-(Q - EQ) / sqrt(varQ))
      }
      pout
    }) 
  }))  
    
  comp.p <- apply(gt@res, 1, function(x) {
    m <- x["test.n"]
    if ((m == 0) | (m == p) )
      out <- NA
    else {
      ps <- sampling.matrix[setsizes == m,]
      out <- mean(ps < (rep(x["p.val"], times=ndraws)) * (1+10^-6))
    }
    out
  })
  if ("comp.p" %in% colnames(gt@res))
    gt@res[,"comp.p"] <- comp.p
  else
    gt@res <- cbind(gt@res, comp.p)
  
  gt
}

#==========================================================


#==========================================================
# The class "gt.barplot" stores the result of a call to
# geneplot or sampleplot
#==========================================================
setClass("gt.barplot", representation(
    res = "matrix", 
    labelsize = "numeric",
    drawlabels = "logical",
    legend = "vector")
    )

#==========================================================
setMethod("result", "gt.barplot",
            function(object) 
{
    res <- object@res
    colouring <- factor(res[,"up"])
    levels(colouring) <- object@legend[4:3]
    res <- data.frame(res[,1:3], z.score(object), colouring)
    colnames(res) <- c("influence", "expected", "sd", "z-score", "info")
    res
})

#==========================================================
# The names method for "gt.barplot" 
#==========================================================
setMethod("names", "gt.barplot", 
            function(x) 
{
  rownames(x@res)
})      

#==========================================================
# A z.score method for "gt.barplot" 
#==========================================================
if( !isGeneric("z.score") )
    setGeneric("z.score", function(x) standardGeneric("z.score"))


setMethod("z.score", "gt.barplot", 
            function(x) 
{
    res <- x@res
    (res[,"influence"] - res[,"Einf"]) / res[,"sd.inf"]
})      


#==========================================================
setMethod("[", "gt.barplot", 
            function(x, i, j,...,drop) 
{
  if (all(i %in% rownames(x@res)) | all(i %in% 1:nrow(x@res))
    | (is.logical(i) & (length(i) == nrow(x@res)))) {
      x@res <- x@res[i, ,drop = FALSE] 
      x
  } else {
    stop("invalid subscript")
  }
})  

#==========================================================
setMethod("length", "gt.barplot", 
            function(x) 
{
  nrow(result(x))
})            

#==========================================================
setMethod("scale", "gt.barplot", 
            function(x, center = FALSE, scale = TRUE) 
{
  if (center) {
    x@res[,"influence"] <- x@res[,"influence"] - x@res[,"Einf"]
    x@res[,"Einf"] <- x@res[,"Einf"] - x@res[,"Einf"]
  }
  if (scale) {
    norm <- x@res[,"sd.inf"]
    norm[norm == 0] <- 1
    x@res[,"influence"] <- x@res[,"influence"] / norm
    x@res[,"Einf"] <- x@res[,"Einf"] / norm
    x@res[,"sd.inf"] <- x@res[,"sd.inf"] / norm
    x
  }
})            

#==========================================================
setMethod("show", "gt.barplot", 
            function(object) 
{
  show(result(object))        
})

#==========================================================
setMethod("plot", "gt.barplot",
  function(x,y,...) {
    plot.gt.barplot <- function(x, genesubset, drawlabels, labelsize, ...) {
      if (!missing(genesubset))
        x <- x[genesubset]
      if (missing(drawlabels))
        drawlabels <- x@drawlabels
      if (missing(labelsize))
        labelsize <- x@labelsize
      influence <- x@res[,"influence"]
      Einf <- x@res[,"Einf"]
      sd.inf <- x@res[,"sd.inf"]
      up <- x@res[,"up"]
      legend <- rownames(x@res)
      m <- length(influence)
      rangebars <- max(0,influence) - min(0,influence)
      minplot <- min(0,influence) 
      maxplot <- max(0,influence) + 0.2 * rangebars
      if (drawlabels & !is.null(legend)){
      # check for space in margin of plot
        plot.new()  
        labwidth <- max(strwidth(legend,"inches",labelsize))
        margins <- par("mai")
        par(new = TRUE, "mai" = c(max(margins[1], labwidth*1.3), margins[2:4]))
        plot( 0, xlim = c(1/2, m+1/2), ylim = c(minplot, maxplot), col = 0, xlab = "", ylab = "influence", xaxt="n",...)
        axis(1, 1:length(legend), legend, cex.axis = x@labelsize, las=2)
        par("mai"=margins)
      }else
        plot( 0, xlim = c(1/2, m+1/2), ylim = c(minplot, maxplot), col = 0, xlab = "nr", ylab = "influence", ...)
      # plot influence bars
      if (m <= 300) {
        rect(xleft = 1:m - 0.4, xright = 1:m + 0.4, ybottom = rep(0,times=m), ytop = influence, col = (up+2), border=0 )
      } else {
        lines(1:m, influence, lwd = 600 / m, type = 'h', col = (up+2))
      }
      # plot Einf reference line
      lines( (1:(m+1))-1/2, c(Einf, Einf[m]), type = "s" )
      # plot sd.inf marks
      nlines <- trunc((influence - Einf) / sd.inf)
      nlines[sd.inf < 10^-3] <- 0
      for (index in 1:max(abs(nlines))) {
        high <- Einf + sign(nlines) * apply(cbind(abs(index), abs(nlines)), 1, min) * sd.inf
        segments( 1:m - 0.4, high, 1:m + 0.4, high)
      }
      # write the legend
      legend(1/2, maxplot, x@legend[1:2], fil = c(3,2))
      invisible(NULL)
    }
    plot.gt.barplot(x, ...)
  }
)


#==========================================================
# Geneplot plots the influence of each gene on the outcome 
#   of the test statistic
# See help(geneplot) for details
#==========================================================

geneplot <- function(gt, geneset, genesubset, scale = FALSE, drawlabels = TRUE, labelsize = 0.6, ...)
{
    # check correct input of gt
    if ( !is(gt, "gt.result"))
        stop("geneplot should be applied to a globaltest result", call. = FALSE)
        
    # extract the right test.genes vector
    if (!missing(geneset))
      gt <- gt[geneset]
    if (length(gt) > 1)
      stop("more than one geneset selected", call. =FALSE)
    if (gt@res[1, "test.n"] == 0)
      stop("empty pathway", call. = FALSE)
    test.genes <- gt@test.genes[[1]]
    if ( !missing(genesubset)) {
      if ( is.character(genesubset) ) {
        test.genes <- intersect(match(genesubset, colnames(gt@X)), test.genes)
        if ( length(test.genes) == 0 ) 
          stop("genesubset is not a subset of the selected geneset")
      }else{
        if (all(genesubset %in% 1:length(test.genes)))
          test.genes <- test.genes[genesubset]
        else
          stop("genesubset is not a subset of the selected geneset")
      }
    }

    # calculate influence per gene and expected influence
    X <- gt@X[,test.genes, drop = FALSE]
    m <- ncol(X)
    Y <- gt@Y
    n <- length(Y)
    XY <- t(X) %*% Y  
    influence <- XY * XY
    up <- (sign(XY) == 1)
    adjX <- t(gt@IminH) %*% X
    if (gt@model == 'logistic') {
      if (gt@df.adjust == 1) {
        tussen <- (adjX * adjX)
        Einf <- colSums(tussen)
        trRR <- Einf * Einf
        tr2R <- Einf * Einf
        trR2 <- colSums(tussen * tussen)
        mu1 <- mean(Y < 0) 
        mumu <- mu1 * (1 - mu1)
        K <- ( 1 - 6 * mu1 + 6 * mu1^2 ) / mumu
        varinf <- K * ( trR2 - tr2R / n ) + 2 * trRR - 2 * tr2R / (n-1)
      } else {
        # in this case Y is not standardized
        Einf <- colSums( sweep(adjX * adjX, 1, gt@pars[["mu2"]], "*") )
        trRR <- Einf * Einf
        mu <- abs(Y)
        mu2 <- gt@pars[["mu2"]]
        mu4 <- mu * (1-mu)^4 + (1-mu) * mu^4
        varinf <- t(adjX * adjX * adjX * adjX) %*% (mu4 - 3 * (mu2 * mu2)) + 2 * trRR
      }
    }
    if (gt@model == 'linear') {
      Einf <- colSums(adjX * adjX)
      trRR <- Einf * Einf
      tr2R <- Einf * Einf
      varinf <- (2 / (n - gt@df.adjust + 2)) * ( (n - gt@df.adjust) * trRR - tr2R )
    }   
    if (gt@model == 'survival') {
      # reconstruct matrices
      times <- gt@pars[["times"]]
      d <- gt@pars[["d"]]
      expci <- gt@pars[["expci"]]
      dtimes <- unique(times[d == 1]) 
      nd <- length(dtimes)
      matrixO <- outer(times, dtimes, "==") * matrix(d, n, nd)
      atrisk <- outer(times, dtimes, ">=")
      hazinc <- as.vector(1 / (expci %*% atrisk))
      matrixP <- outer(expci, hazinc, "*") * atrisk
      matrixPO <- matrixP %*% t(matrixO)
      matrixM <- diag(d) %*% outer(times, dtimes, "<") - matrixPO %*% outer(times, dtimes, "<")
      matrixW <- diag(rowSums(matrixPO)) - matrixPO %*% t(matrixPO)
      Einf <- colSums(adjX * (matrixW %*% adjX))
      varinf <- sapply( 1:nd, function(j) {
        pj <- matrixP[,j]
        mj <- matrixM[,j]
        oj <- matrixO[,j]
        part1 <- (diag(n) - outer(rep(1,n), pj)) %*% (adjX * adjX)
        part2 <- sweep((diag(n) - outer(rep(1,n), pj)) %*% adjX, 2, as.vector(t(adjX) %*% (mj - pj)), "*")
        tj <- part1 + 2 * part2
        as.vector(pj %*% (tj * tj)) * sum(oj)
      })
      if (is.matrix(varinf))
        varinf <- rowSums(varinf)
      else
        varinf <- sum(varinf)
      tiejs <- (1:nd)[colSums(matrixO) > 1]
      if (length(tiejs) > 0) {
        influence <- influence - rowSums( sapply(tiejs, function(j) {
          pj <- matrixP[,j]
          oj <- matrixO[,j]
          matrixtie <- outer(oj, oj) - diag(oj)
          pjX <- as.vector(pj %*% X)
          colSums(X * (matrixtie %*% X)) - 2 * pjX * as.vector(rowSums(matrixtie) %*% X) + sum(matrixtie) * pjX * pjX 
        }))
      }
      # reconstruct expected variance and use to standardize
      R <- adjX %*% t(adjX) / m
      tussen <- matrix(diag(R), n, nd) + 2 * R %*% (matrixM - matrixP)
      matrixT <- tussen - matrix(colSums(matrixP * tussen), n, nd, byrow = TRUE)
      varQ <- sum(matrixPO * (matrixT * matrixT) %*% t(matrixO) )
      norm <- sqrt(varQ)
      influence <- (influence - Einf) / norm
      Einf <- 0
      varinf <- varinf / (norm * norm)
    }
    sd.inf <- sqrt(varinf)

    # Output: gt.barplot object
    res <- cbind(influence, Einf, sd.inf, up)
    colnames(res) <- c("influence", "Einf", "sd.inf", "up")
    rownames(res) <- colnames(X)
    if (gt@model == 'linear')
      legend <- c("positive correlation with residual Y", "negative correlation with residual Y", "+", "-")
    if (gt@model == 'logistic')
      legend <- c(paste("higher expression in", gt@levels[2], "samples"), 
        paste("higher expression in", gt@levels[1], "samples"), paste("high in", gt@levels[2]), 
        paste("high in", gt@levels[1]) )
    if (gt@model == 'survival')
      legend <- c(paste("positively associated with survival"), 
        paste("negatively associated with survival"), "+", "-")
    
    gtbar <- new("gt.barplot",     
      res = res,
      labelsize = labelsize, 
      drawlabels = drawlabels,
      legend = legend)
    if (scale)
      gtbar <- scale(gtbar)
    plot(gtbar)
    invisible(gtbar)
}
#==========================================================


#==========================================================
# Sampleplot plots the influence of each sample on the outcome 
#   of the test statistic
# See help(sampleplot) for details
#==========================================================
sampleplot <- function(gt, geneset, samplesubset, scale = FALSE, drawlabels = TRUE, labelsize = 0.6,...)
{
    # check correct input of gt
    if ( !is(gt, "gt.result"))
        stop("geneplot should be applied to a globaltest result", call. = FALSE)
        
    # extract the right test.genes vector
    if (!missing(geneset))
      gt <- gt[geneset]
    if (length(gt) > 1)
      stop("more than one geneset selected", call. =FALSE)
    if (gt@res[1, "test.n"] == 0)
      stop("empty pathway", call. = FALSE)
    test.genes <- gt@test.genes[[1]]

    # reconstruct X and Y
    X <- as.matrix(gt@X[,test.genes])
    n <- nrow(X)
    if (missing(samplesubset))
      samplesubset <- 1:n
    else {  
      if (all(samplesubset %in% rownames(X)))
        samplesubset <- match(samplesubset, rownames(X))
      else
        if (!all(samplesubset %in% 1:n) & !all(samplesubset %in% rownames(X)))
          stop("samplesubset should contain names or numbers of samples", call. = FALSE)
    }
    m <- ncol(X)
    Y <- gt@Y
    
    # calculate influence per sample and expected influence
    if (gt@model != 'survival') {
      XXY <- t(gt@IminH) %*% X %*% t(X) %*% Y
      influence <- (n/m) * Y * XXY
      up <- (sign(Y) == 1) 
      if ((gt@model == 'logistic') & (gt@df.adjust > 1)) {
        mu2 <- abs(Y) * (1 - abs(Y))
        adjX <- t(gt@IminH) %*% X
        R <- adjX %*% t(adjX)
        RR <- (R * R) %*% diag(mu2)
      } else {
        adjX <- t(gt@IminH) %*% X
        R <- adjX %*% t(adjX)
        RR <- R * R
      }
      Einf <- (n/m) * Y * Y * rowSums(adjX*adjX)
      RR <- matrix(t(RR)[diag(n) == 0], n, n-1, byrow = TRUE)
      varinf <- (n/m)^2 * Y * Y * rowSums(RR)
      sd.inf <- sqrt(varinf)
    } else {
      # survival model
      adjX <- t(gt@IminH) %*% X
      R <- adjX %*% t(adjX) 
      times <- gt@pars[["times"]]
      d <- gt@pars[["d"]]
      expci <- gt@pars[["expci"]]
      dtimes <- unique(times[d == 1]) 
      nd <- length(dtimes)
      matrixO <- outer(times, dtimes, "==") * matrix(d, n, nd)
      atrisk <- outer(times, dtimes, ">=")
      hazinc <- as.vector(1 / (expci %*% atrisk))
      matrixP <- outer(expci, hazinc, "*") * atrisk
      matrixPO <- matrixP %*% t(matrixO)
      matrixM <- diag(d) %*% outer(times, dtimes, "<") - matrixPO %*% outer(times, dtimes, "<")
      matrixW <- diag(rowSums(matrixPO)) - matrixPO %*% t(matrixPO)
      influence <- rowSums((adjX %*% t(X)) * (outer(Y,Y) - matrixW))
      up <- (sign(Y) == 1) 
      eye <- diag(n)
      Evarinf <- sapply(1:n, function(i) {
        out <- numeric(2)
        Ri <- R[i,]
        mPi <- matrixP[i,]
        mMi <- matrixM[i,]
        part1 <- matrix(mMi - mPi, n, nd, byrow = TRUE) * outer(Ri, as.vector(Ri %*% matrixP), "-")
        part2 <- matrix(Ri %*% (matrixM - matrixP) + rep(Ri[i],nd), n, nd, byrow=TRUE) * (matrix(eye[i,],n,nd) - matrix(mPi,n,nd, byrow = TRUE))
        matrixK <- part1 + part2
        minusi <- !as.logical(eye[i,])
        di <- as.logical(matrixO[i,])
        out[1] <- ifelse(any(di), matrixK[i,di], 0) + sum(matrixK[minusi,] * matrixP[minusi,])
        matrixK[i,] <- 0
        out[2] <- sum((matrixP %*% t(matrixO)) * ((matrixK * matrixK) %*% t(matrixO)))
        out
      })
      Einf <- Evarinf[1,]
      sd.inf <- sqrt(Evarinf[2,])
      tiecorrect <- sum( sapply(1:nd, 
        function(j) {
          oj <- matrixO[,j]
          if (sum(oj) == 1)
            rep(0, n)
          else { 
            XX <- X %*% t(X) / n
            matrixtie <- outer(oj, oj) - diag(oj)
            pj <- matrixP[,j]
            pjXX <- as.vector(XX %*% pj)
            tc <- diag(matrixtie %*% XX) - rowSums(matrixtie) * pjXX - pj * rowSums(matrixtie) %*% XX + pj * pjXX * sum(matrixtie)
            tc
          }
        }
      ))
      influence <- influence - tiecorrect
      # reconstruct expected variance and use to standardize
      tussen <- matrix(diag(R), n, nd) + 2 * R %*% (matrixM - matrixP)
      matrixT <- tussen - matrix(colSums(matrixP * tussen), n, ncol(matrixP), byrow = TRUE)
      varQ <- sum((matrixP %*% t(matrixO)) * ((matrixT * matrixT) %*% t(matrixO)))
      norm <- sqrt(varQ) / n
      influence <- influence / norm
      Einf <- Einf / norm
      sd.inf <- sd.inf / norm
    }
    
    # Output: gt.barplot object
    res <- cbind(influence - Einf, 0, sd.inf, up)
    colnames(res) <- c("influence", "Einf", "sd.inf", "up")
    rownames(res) <- rownames(X)
    if (gt@model == 'linear')
      legend <- c("positive residual Y", "negative residual Y", "+", "-")
    if (gt@model == 'logistic')
      legend <- c(paste(gt@levels[2], "samples"), paste(gt@levels[1], "samples"), gt@levels[2], gt@levels[1])
    if (gt@model == 'survival')
      legend <- c("late event time or censored", "early event time", "late", "early")
    
    gtbar <- new("gt.barplot",     
      res = res,
      labelsize = labelsize, 
      drawlabels = drawlabels,
      legend = legend)
    gtbar <- gtbar[samplesubset]
    if (scale)
      gtbar <- scale(gtbar)
    plot(gtbar)
    invisible(gtbar)
}
#==========================================================

#==========================================================
# The function regressionplot allows the evaluation of 
#   possibly outlying samples.
# See help(regressionplot) for details
#==========================================================

regressionplot <- function(gt, geneset, sampleid, ...)
{

    # check correct input of gt
    if ( !is(gt, "gt.result"))
        stop("geneplot should be applied to a globaltest result", call. = FALSE)
        
    # extract the right test.genes vector
    if (!missing(geneset))
      gt <- gt[geneset]
    if (length(gt) > 1)
      stop("more than one geneset selected", call. =FALSE)
    if (gt@res[1, "test.n"] == 0)
      stop("empty pathway", call. = FALSE)
    test.genes <- gt@test.genes[[1]]

    # recreate matrix R and find S = Y %o% Y
    X <- as.matrix(gt@X[,test.genes])
    Y <- gt@Y
    R <- as.vector(X %*% t(X))
    S <- as.vector(Y %o% Y)
    n <- length(Y)
    
    # Check correct input of sampleid
    if ( missing(sampleid) ) {
      sampleid <- NULL
    }else{
      if (!( all(sampleid %in% 1:n) | all(sampleid %in% rownames(X)) ))
        stop("Option sampleid incorrect", call. = FALSE)
    }
            
    # Extract relevant entries from S and R
    samples <- rep(FALSE,times = n)
    names(samples) <- rownames(X)
    samples[sampleid] <- TRUE
    samples <- outer(samples, samples, "|")
    lowertriangle <- outer(1:n, 1:n, ">")
    selection <- as.vector(lowertriangle)
    subselection <- as.vector( lowertriangle & samples )
    Rall <- R[selection]
    Sall <- S[selection]
    Rsub <- R[subselection]
    Ssub <- S[subselection]
    Rrest <- R[selection & !subselection]
    Srest <- S[selection & !subselection]
    
    # Draw the plots
    plot(Sall, Rall, xlab = "Covariance between outcomes", ylab = "Covariance between expression profiles", col = 0,...)
    if (length(Rrest) > 0){
      points(Srest, Rrest,col = 4, cex=0.5)
      abline(lm(Rall ~ Sall), col = 4)
    }
    if (length(Rsub) > 0){
      points(Ssub, Rsub, col = 2, pch = 4, cex = 2)
      abline(lm(Rsub ~ Ssub), col = 2, lty = 2)
    }
    
    # No output
    invisible(NULL)
}   
#==========================================================


#==========================================================
# Function Checkerboard visualizes a result of globaltest
#   by visualizing the covariance matrix between the different samples
# See help(checkerboard) for details
#==========================================================

checkerboard <- function(gt, geneset, sort = TRUE, drawlabels = TRUE, labelsize = 0.6,...)
{   
    # check correct input of gt
    if ( !is(gt, "gt.result"))
        stop("geneplot should be applied to a globaltest result", call. = FALSE)
        
    # extract the right test.genes vector
    if (!missing(geneset))
      gt <- gt[geneset]
    if (length(gt) > 1)
      stop("more than one geneset selected", call. =FALSE)
    if (gt@res[1, "test.n"] == 0)
      stop("empty pathway", call. = FALSE)
    test.genes <- gt@test.genes[[1]]
    if (!is.logical(sort))
        stop("The option 'sort' should be either TRUE or FALSE", call. = FALSE)

    # recreate matrix R
    X <- as.matrix(gt@X[,test.genes])
    Y <- gt@Y
    R <- X %*% t(X)
        
    # Sort Y if needed and return new samplenrs
    n <- length(Y)
    perm <- sort.list(Y)
    rperm <- sort.list(Y[n:1])
    if ( any(perm != 1:n ) & any(rperm != 1:n) & (sort)){
        newsamplenrs <- matrix( c(1:n, sort.list( perm )), n, 2 ) 
        label <- "sorted samplenr"
        R <- R[perm,perm]
    }else{
        sort = FALSE
        label <- "samplenr"
        newsamplenrs <- matrix( c(1:n, 1:n), n, 2 ) 
    }    
    colnames(newsamplenrs) <- c("samplnr.old", "samplenr.new")
    rownames(newsamplenrs) <- rownames(X)

    # Calculate median non-diagonal element of R
    lowertriangle <- outer( 1:n, 1:n, ">" )
    med <- median(R[lowertriangle])
    # Draw plot
    if (is.null(labelsize))
      labelsize<-par("cex.axis")
    if (drawlabels & !is.null(rownames(newsamplenrs))){
      legend<-rownames(newsamplenrs)[sort(newsamplenrs[,2],index.return=TRUE)$ix]   
      # check for space in margin of plot
      plot.new()  
      labwidth<-max(strwidth(legend,"inches",labelsize))
      margins<-par("mai")
      par(new=TRUE,"mai"=c(max(margins[1],labwidth*1.3),max(margins[2],labwidth*1.3), margins[3:4]))
      image(x = 1:n, y = 1:n, z = R>med, col = rainbow(2, v = c(0,1), s = c(1,0) ), ylab = "", xlab = "", 
          yaxt="n", xaxt = "n", ...)
      axis(2,1:length(legend), legend, cex.axis=labelsize, las=2)
      axis(1,1:length(legend), legend, cex.axis=labelsize, las=2)
      par("mai"=margins)
    }else
      image(x = 1:n, y = 1:n, z = R>med, col = rainbow(2, v = c(0,1), s = c(1,0) ), xlab = label, ylab = label, 
          lab = c(n,n,50/n), ...)
    par(pty = "s")
    invisible(newsamplenrs)
}
#==========================================================

#==========================================================
# .First.lib is called when the package is loaded
#==========================================================

.First.lib <- function(libname, pkgname, where)
{
    if (missing(where)) {
        where <- match(paste("package:", pkgname, sep=""), search())
        if(is.na(where)) {
            warning(paste("Not a package name: ",pkgname))
            return()
        }
        where <- pos.to.env(where)
    }
  invisible(NULL)
}
#==========================================================
