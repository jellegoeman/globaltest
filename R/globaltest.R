#==========================================================                                        
# Function "globaltest" performs the global test on (a list of) 
#    subsets of the data
# and stores the result in a gt.result object
#==========================================================

globaltest <- function(X, Y, genesets,
                        model, levels,
                        d, event = 1,
                        adjust, 
                        method = c("auto", "asymptotic", "permutations", "gamma"),
                        nperm = 10^4,
                        scaleX = TRUE, 
                        accuracy = 50, ...)

{

  # We coerce the input X, Y, d, model, adjust, levels, into the format:
  # matrix eX
  # data.frame pData
  # formula ff

  # 0: check whether the input is in the permitted classes
  if (!(class(X) %in% c("matrix","ExpressionSet"))) 
    stop("X should be of class ExpressionSet or matrix", call. = FALSE)
  if (!((class(Y) %in% c("factor", "formula")) || is.vector(Y) ))
    stop("Y should be of class formula, vector, factor or character", call. = FALSE)
  if (!(missing(genesets) || is.vector(genesets)))
    stop("genesets should be of class vector or list", call. = FALSE)
  if (!(missing(levels) || is.vector(levels)))
    stop("levels should be of class vector", call. = FALSE)
  if (!(missing(d) || is.vector(d))) 
    stop("d should be of class vector", call. = FALSE)
  if (length(event) > 1)
    stop("event should contain a single value", call. = FALSE)
  if (!(missing(adjust) || (class(adjust) %in% c("character","formula", "vector", "data.frame")))) 
    stop("adjust should be of class character, formula, vector or data.frame", call. = FALSE)
  if (!(class(scaleX) %in% c("logical"))) 
    stop("scaleX should be either TRUE or FALSE", call. = FALSE)
  if (accuracy <= 1)
    stop("accuracy should be greater than 1")
 
  
  # 1: extract the expression matrix:
  if (is(X, "ExpressionSet")) {
    eX <- exprs(X)
  } else {
    if (is.data.frame(X) || is.vector(X)) {
      X <- as.matrix(X)
    }
    if (is.matrix(X)) {
      eX <- X
      # maybe X has to be transposed later
    }
  }
  
  # 2: extract the phenoData from X, adjust and/or Y
  pDataNamesSupplied <- TRUE
  pData <- NULL
  if (is(X, "ExpressionSet")) {
    pData <- pData(X)
  }
  if (!missing(adjust) && is.data.frame(adjust)) {
    if (is.null(pData)) {
      pData <- adjust
    } else {
      if (length(intersect(names(pData), names(adjust))) > 0) {
        warning("Duplicate variable names in pData(X) and adjust", call. = FALSE)
      }
      if ((nrow(adjust) != nrow(pData)) || !all(rownames(pData) == rownames(adjust))) {
        warning("sample names in pData(X) and adjust do not match", call. = FALSE)
      }
      pData <- data.frame(pData, adjust)
      adjust <- names(adjust)
    }
  }
  if ((!is(Y, "formula")) && (length(Y) > 1)) {
    if (is.null(pData)) {
      pDataNamesSupplied <- !is.null(names(Y))
      pData <- data.frame(Y)
    } else {
      if ("Y" %in% names(pData)) {
        stop("Variable name Y not allowed in data.frames pData(X) or adjust", call. = FALSE)
      }
      if (length(Y) != nrow(pData)) {
        stop("The number of samples in Y does not match that in pData(X) and/or adjust", call. = FALSE)
      }
      if (!is.null(names(Y)) && !all(names(Y) == rownames(pData))) {
        warning("The sample names in Y do not match those in pData(X) and/or adjust", call. = FALSE)
      }
      pData <- data.frame(pData, Y)
    }
  }
  if (!missing(d) && is.vector(d) && (length(d) > 1)) {
    if ("d" %in% names(pData)) {
      stop("Variable name d not allowed in data.frames pData(X) or adjust", call. = FALSE)
    }
    if (length(d) != nrow(pData)) {
      stop("The number of samples in d does not match that in Y, pData(X) and/or adjust", call. = FALSE)
    }
    if (!is.null(names(d)) && !all(names(d) == rownames(pData))) {
      warning("The sample names in d do not match those in Y, pData(X) and/or adjust", call. = FALSE)
    }
    pData = data.frame(pData, d)
  }
    
  # 3: check for sample size conflict between eX and pData
  if ((!is.null(pData)) && (ncol(eX) != nrow(pData))) {
    if (nrow(eX) == nrow(pData)) {
      eX <- t(eX)
    } else {
      stop("Number of samples in X or exprs(X) conflicts with number of samples in Y, pData(X) or adjust", call. = FALSE)
    }
  }
  n <- ncol(eX)
  p <- nrow(eX)
  if (is.null(pData)) {
    pData <- as.data.frame(matrix(,n,0))
    pDataNamesSupplied <- FALSE
  } 

  # 4: Synchronize sample names between eX and pData
  if (pDataNamesSupplied && !is.null(colnames(eX))) {
    if (!all(rownames(pData) == colnames(eX))) {
      warning("Sample names in X or exprs(X) conflict with those in Y, pData(X) or adjust", call. = FALSE)
    }
  } else {
    if (is.null(colnames(eX))) {
      colnames(eX) <- rownames(pData)
    } else {
      if (!pDataNamesSupplied) {
        rownames(pData) <- colnames(eX)
      }
    }
  }
  
  # 5: remove the input X to free memory
  remove(X)
  
  # 6: Remove redundant levels from factor variables in pData
  if (ncol(pData) > 0) {
    pData <- data.frame(lapply(as.list(pData), function(x) { if (is.factor(x)) factor(x) else x}))
  }
  rownames(pData) <- colnames(eX)

  # 7: Extract the formula object
  if (is(Y, "formula")) {
    ff <- Y
    remove(Y)
  } else {
    if (!missing(adjust) && is(adjust, "formula")) {
      ff <- adjust
    } else {
      if (is.character(Y) && length(Y) == 1) {
        ffstart = paste(Y, "~")
      } else {
        ffstart = "Y ~"
      }
      if (!missing(adjust)) {
        if (is.character(adjust)) {
          ff <- formula(paste(ffstart, paste(adjust, collapse = "+")))
        } else if (is.data.frame(adjust)) {
          ff <- formula(paste(ffstart, paste(names(adjust), collapse = "+")))
        }
      } else {
        ff <- formula(paste(ffstart, "1"))
      }
    }
  }
  absent <- all.vars(ff)[sapply(all.vars(ff), function(vr) {
    (length(find(vr)) == 0) && !(vr %in% names(pData))
  })]
  if (length(absent) > 0) {
    stop(paste("Variable(s)", paste(absent, collapse = ", "), "not found"), call. = FALSE)
  }

  # 8: check for missing values in pData
  outcome <- eval(ff[[2]], envir = pData)
  if (any(sapply(eval(attr(terms(ff), "variables"), envir = pData), function(x) any(is.na(x)))) || any(is.na(outcome)))
      stop("Missing values only allowed in the gene expressions.", call. = FALSE)
  
  
  # 9: determine the model from the input if missing
  if (!missing(model)) {
    model <- match.arg(model, c('linear','logistic','survival'))
  }
  if (missing(model)) {
    if (is(outcome, "factor") || !missing(levels) || (length(unique(outcome)) <= 2)) {
      model <- "logistic"
    } else if (is(outcome, "Surv") || !missing(d)) {
      model <- "survival"
    } else if (is(outcome, "numeric")) {
      model <- "linear"
    } else {
      stop("The model could not be determined from the input. Please specify the model.", call. = FALSE)
    }
  }
  
  
  # 10: adept the formula in case of a survival or logistic model
  if ((model == "survival") && !is(outcome, "Surv")) {
    if (!missing(d)) {
      if (length(d) == 1) {
        ff <- formula(paste("Surv(", ff[[2]], ",", d, " == event)~", ff[[3]])) 
      } else {
        ff <- formula(paste("Surv(", ff[[2]], ", d == event)~", ff[[3]]))
      }
    } else { 
        ff <- formula(paste("Surv(", ff[[2]], ")~", ff[[3]]))
    } 
  }  
  if ((model == "logistic") && is.numeric(outcome) && !all(outcome %in% 0:1)) {
    ff <- formula(paste("factor(", ff[[2]], ")~", ff[[3]]))
  }
     
  # 11: Preparation of eX and pData in the logistic model using option levels
  if (model == 'logistic') {
    if (missing(levels)) {
      levels <- NULL
      levels <- levels(factor(outcome))
      if (length(levels) == 1) {
        stop("Y has only one value.", .call = FALSE)
      }
    } else {
      if (!all(levels %in% levels(factor(outcome)))) {
        absent <- levels[!(levels %in% levels(factor(outcome)))]
        stop(paste("Levels", paste(absent, collapse = ", "), "not present in Y"), call. = FALSE)
      }
    }
    if (length(levels) == 1) {
      # in this case we compare one level against the others
      # a new variable must be created
      newnameY <- paste(as.character(ff[[2]]), ".", make.names(as.character(levels[1])), sep = "", collapse = "")
      newY <- rep(levels[1], n)
      newY[outcome != levels[1]] <- "other"
      newY <- factor(newY, levels = c(levels[1], "other"))
      pData <- data.frame(pData, newY)
      names(pData)[length(pData)] <- newnameY
      ff <- formula(paste(newnameY, "~", as.character(ff[[3]])))
    } else if (length(levels) < length(levels(factor(outcome)))) { 
      eX <- eX[,outcome %in% levels,drop = FALSE]
      vars <- all.vars(ff)[!(all.vars(ff) %in% names(pData))]
      if (length(vars) > 0) {
        add2pData <- as.data.frame(sapply(vars, function(var) {
          eval(as.name(var))
        }))
        pData <- data.frame(pData, add2pData)
      }
      pData <- pData[outcome %in% levels,, drop = FALSE]
      n <- ncol(eX)
    }
    if (length(levels) > 2) {
      model <- "multinomial"  
      g <- length(levels)
    }
  } 
                                   
  # 12: Fit the adjustmodel:
  if (model == "logistic") {
    fit <- glm(ff, data = pData, family = binomial, x = TRUE, y = TRUE, na.action = 'na.fail', 
      control = glm.control(maxit = 20))
  } else if (model == "linear") {
    fit <- lm(ff, data = pData, x = TRUE, y = TRUE, na.action = 'na.fail')
  } else if (model == "survival") {
    fit <- coxph(ff, data = pData, x = TRUE, y = TRUE, na.action = 'na.fail')
  } else if (model == "multinomial") {
    fit <- mlogit(ff, data = pData, control = glm.control(maxit = 20))
  }
  if (!is.null(coefficients(fit)) && any(is.na(coefficients(fit))))
    stop("The null model is singular.", call. = FALSE)
  if (!is.null(iterations(fit)) && iterations(fit) == 20)
    warning("The null model failed to converge.", call. = FALSE)
  if ((model == "survival") && (all(residuals(fit) == 0))) {
    stop("There are no events.", call. = FALSE)
  }


  # 13: Calculate IminH
  IminH <- .getIminH(fit)
  adjusted <- !is.null(IminH)


  # 14: Impute missing values of X and rescale X       
  row.meanX <- rowMeans(eX, na.rm = TRUE)
  eX <- eX - row.meanX %o% rep(1, times=n)
  eX[is.na(eX)] <- 0
  #if (correlation) {
  #  eX <- eX / apply(eX, 1, sd) %o% rep(1, times=n)
  #}
  if (scaleX && !(model %in% c("survival", "multinomial"))) {
    if (adjusted) {
      adjX <- eX %*% IminH
    } else {
      adjX <- eX
    }
    if ((model == "logistic") && adjusted) {
      norm <- sqrt(sum((adjX * adjX) %*% fit$weights) / (10 * p))
    } else {
      norm <- sqrt(sum(adjX * adjX) / (10 * p))
    }
    remove(adjX)
    eX <- eX / norm
  }
  if (scaleX && (model == "multinomial")) {
    if (adjusted) {
      norm <- sqrt(.EQall(eX, IminH, fit) / 10)
    } else {
      mu <- fitted.values(fit)[1,]
      matrixW <- -mu %o% mu + diag(mu)
      XX <- crossprod(eX) / p
      norm <- sqrt(sum(diag(matrixW)) * sum(diag(XX)) / 10)
    }
    eX <- eX / norm
  }



  # 15: Coerce genesets into a list of vectors of numbers and check correct input
  if ( missing(genesets) ) {
    if ("test.genes" %in% names(list(...))) {
      genesets <- list(...)[["test.genes"]]
    } else {
      genesets <- list(all = 1:p)
    }
  }
  
  if ( !is.list(genesets) ) 
    genesets <- list(genesets)

  genesets <- lapply(genesets, function(tg) { 
    if ( !is.vector(tg) & !is.null(tg) )
      stop("genesets should be a (list of) vector(s)", call. = FALSE)

    if ( (length(tg) == p) && all(tg %in% c(0,1)) ) {
      test.names <- names(tg)[tg == 1]
      if (  !is.null(rownames(eX) )  &  !is.null(test.names) ) {
        if( any( rownames(eX)[tg] != test.names ) ) 
          warning("Gene names in X inconsistent with gene names in genesets.", call. = FALSE)
      }
      tg <- (1:p)[tg == 1]
    } else {
      if (is.character(tg))  {
        tg <- match(tg, rownames(eX))
      } else {  
        if ( all(tg %in% 1:p) ) {
          test.names <- names(tg)
          names(tg) <- NULL
          if (  !is.null(rownames(eX) )  &  !is.null(test.names) ) {
            if( any( rownames(eX)[tg] != test.names ) ) 
            warning("Gene names in X inconsistent with gene names in genesets.", call. = FALSE)
          }
        } else {
          stop("option 'genesets' should be a (list of) vector(s) of gene names or numbers.", call. = FALSE)
        }
      } 
    }  
    tg
  })
  path.n <- sapply(genesets, length) 
  
  # 16: Throw away all genes in genesets that are not on the chip
  genesets <- lapply(genesets, function(tg) if (!is.null(tg)) unique(tg[!is.na(tg)]) else tg ) 
  test.n <- sapply(genesets, length)


  # 17: Prepare the return object of type "gt.result"
  if (!adjusted)
    IminH <- matrix(,0,0)
  nullres <- cbind(Genes = path.n, Tested = test.n)
  rownames(nullres) <- names(genesets)
  Perms <- matrix(,length(genesets),0)
  rownames(Perms) <- names(genesets)
  gt <- new("gt.result",
    res = nullres,
    eX = eX,
    genesets = genesets,
    fit = list(fit),
    IminH = IminH,
    PermQs = Perms
  )

  # 18: Check input of "method"
  method <- match.arg(method)
  if (method == "auto") {
    if (!((model %in% c("linear", "logistic")) && adjusted) && (.nPerms(gt) < nperm)) {
      method <- "permutations"
    } else {
      method <- "asymptotic"
    }
  }
  if ((model == "survival") && (method == "gamma")) {
    warning("No gamma approximation for a survival model. Asymptotic normality used.", call. = FALSE)
    method <- "asymptotic"
  }
  if (method == "gamma") {
    warning("The gamma method has been deprecated. It will be removed in a future version.")
  }
  

  # 19: Calculate the test results for each geneset and add them to gt
  if ((method == "gamma") || (method == "permutations")) {
    gt@method <- 2
  } else {
    gt@method <- 3
  }
  res <- .globaltest(gt, accuracy = accuracy)
  if (ncol(res) > 0) {
    colnames(res) <- c("Statistic Q","Expected Q","sd of Q","P-value")
    gt@res <- cbind(gt@res, res)
  }
  if (method == "permutations") {
    gt <- permutations(gt, nperm = nperm)
  }

  # 20: return
  gt
}
#==========================================================

  
  
