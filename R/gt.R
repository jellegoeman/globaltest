gt <- function(response, alternative, null, data, test.value,    
            model = c("linear", "logistic", "cox", "poisson", "multinomial"), levels,
            directional = FALSE, standardize = FALSE, permutations = 0, 
            subsets, weights, alias, x = FALSE, trace) {
  
  # store the call
  call <- match.call()
  
  # avoid conflict between "levels" input and "levels" function
  if (missing(levels)) levels <- NULL
                       
  # data default
  if (missing(data)) {
    if ((!missing(alternative)) && (is(alternative, "ExpressionSet"))) {
      data <- pData(alternative)
    } else
      data <- NULL
  }
  if (is.matrix(data))  
    data <- data.frame(data)

  # evaluate response, which may be one of the colnames of data
  response <- eval(call$response, data, parent.frame())

  # settle null, alternative and response if response is a formula
  if (missing(null) || is.null(null)) {
    if ((!missing(alternative)) && is(response, "formula"))
      null <- response
    else
      null <- ~1
  } 
  if (missing(alternative))
    if (is(response, "formula"))
      alternative <- response
    else
      stop("argument \"alternative\" is missing, with no default")  
  if (is(response, "formula")) {
    name.response <-  deparse(eval(response)[[2]])
    response <- eval(attr(terms(response, data=data), "variables"), data, environment(response))[[attr(terms(response, data=data), "response")]]
  } else {
    name.response <- deparse(call$response)
  }
  
  # remove redundant levels from factor response
  # and coerce response to factor in case of levels input
  if (is.factor(response) || !is.null(levels))
    response <- factor(response) 
                                     
  # get the model
  if (missing(model)) {
    if(is(response, "Surv")) 
      model <- "cox"
    else if ((is.factor(response) && length(levels(response)) <= 2) || is.logical(response)) 
      model <- "logistic"
    else if (is.factor(response) && length(levels(response)) > 2)
      model <- "multinomial"
  }
  model <- match.arg(tolower(model), c("linear", "logistic", "cox", "poisson", "multinomial"))

  # if multinomial, coerce to factor and remove redundant level
  if (model=="multinomial" && !is.factor(response))
      response <- factor(response)
  
  # if multinomial, coerce to logistic if possible
  if (model=="multinomial" && nlevels(response)==2)
    model <- "logistic"
  
  # find the sample size
  if (model == "cox") {
    if (attr(response, "type") %in% c("right", "counting"))
      n <- nrow(response)
    else
      stop("survival data of type", attr(response, "type"), "not supported")
  } else {
    n <- length(response)
  }
                
  # remove terms from alternative that are also in null
  if (is(null, "formula") && is(alternative, "formula") && 
        identical(environment(null), environment(alternative))) {
    dup <- attr(terms(alternative, data=data), "term.labels") %in% attr(terms(null, data=data), "term.labels")
    if (all(dup)) stop("all covariates in alternative also in null")  
    if (any(dup)) 
      alternative <- formula(terms(alternative,data=data)[!dup])
  }

  # get null and alternative
  null <- .getNull(null, data, n, model)
  alternative <- .getAlternative(alternative, data, n)

  # Adjust input due to levels argument
  if ((!is.null(levels)) && is.factor(response)) {
    if (!all(levels %in% levels(response)))
      stop("argument \"levels\" does not match levels(response)")
    if (length(levels) > 1) {
      select <- response %in% levels
      response <- factor(response[select], levels=levels)
      alternative <- alternative[select,,drop=FALSE]
      null$null <- null$null[select,, drop=FALSE]
      if (!is.null(null$offset)) 
        null$offset <- null$offset[select]
      if (length(levels) == 2)
        model <- "logistic"
    } else {
      response <- factor(response == levels)
      levels(response) <- c("other", levels)
      model <- "logistic"
    }
  }

  # prepare legends for plots later
  legend <- list()
  if (model == "logistic") {
    if (is.factor(response)) {
      legend$cov <- paste("assoc. with", name.response, "=", rev(levels(response)))
      legend$subj <- paste(name.response, "=", rev(levels(response)))
    } else {
      legend$cov <- paste("assoc. with", name.response, "=", 1:0)
      legend$subj <- paste(name.response, "=", 1:0)
    }
  } else if (model == "multinomial") {
    legend$cov <- paste("assoc. with", name.response, "=", levels(response))
    legend$subj <- paste(name.response, "=", levels(response))
  } else if (model %in% c("linear", "poisson")) {
    legend$cov <- paste(c("pos. assoc. with", "neg. assoc. with"), name.response)
    legend$subj <- paste(c("pos. residual", "neg. residual"), name.response)
  } else { # model = "survival" 
    legend$cov <- c("pos. assoc. with survival", "neg. assoc. with survival")
    legend$subj <- c("late event or censored", "early event")
  }

  #subsets and weigths
  if (!missing(subsets) && !is.list(subsets))
    subsets <- list(subsets)
  many.subsets <- !missing(subsets)
  one.weight <- (!missing(weights)) && (!is.list(weights)) && (length(weights)==ncol(alternative)) && many.subsets
  many.weights <- (!missing(weights)) && (!one.weight)
  if (many.weights && !is.list(weights))
    weights <- list(weights)

  # if a test.value was specified, adjust the offset term
  # NOTE: test.value is interpreted relative to alternative BEFORE standardizaton and weighting
  if (!missing(test.value)) {
    if (length(test.value) == ncol(alternative))
      test.offset <- drop(alternative %*% test.value)
    else if ((length(subsets)==1) && (length(test.value) == length(subsets[[1]])))
      test.offset <- drop(alternative[,subsets[[1]], drop=FALSE] %*% test.value)
    else
      stop("the length of \"test.value\" (",length(test.value), ") does not match the the number of covariates in \"alternative\" (", ncol(alternative), ")")
    if (!is.null(null$offset))
      offset <- null$offset + test.offset
    else      
      offset <- test.offset
  } else
    offset <- null$offset   
  null <- null$null
  if (model == "multinomial" && !is.null(offset))
    stop("offset term and test.value not yet implemented for the multinomial model")

  # conservatively impute missing values in alternative
  all.na <- apply(is.na(alternative), 2, all)
  some.na <- apply(is.na(alternative), 2, any) & !all.na
  if (ncol(null) == 0) {
    if (model == "cox") {
      alternative[,some.na] <- apply(alternative[,some.na, drop=FALSE], 2, function(cov) {
        cov[is.na(cov)] <- mean(cov, na.rm=TRUE)
        cov
      })
    } else
      alternative[is.na(alternative)] <- 0
  } else {
    alternative[,some.na] <- apply(alternative[,some.na, drop=FALSE], 2, function(cov) {
      fit <- lm(cov ~ 0 + null, x = TRUE)
      coefs <- coef(fit)
      coefs[is.na(coefs)] <- 0
      cov[is.na(cov)] <- drop(null %*% coefs)[is.na(cov)]
      cov
    })
    alternative[,all.na] <- 0 
  }
                            
  # deafult alias
  if (missing(alias)) alias <- NULL

  # check weights and subsets lengths
  if (many.subsets && many.weights) {
    if (length(subsets) != length(weights))
      stop("lengths of \"subsets\" and \"weights\" do not match")
    if ((!is.null(names(subsets))) && (!is.null(names(weights))) && (!all(names(subsets)==names(weights))))
      if (is.null(alias))
        alias <- names(weights)
      else
        warning("names of subsets and weights do not match")
  }
  
  # make sure subsets is a list of names, compatible with colnames(alternative)
  if (many.subsets && gt.options()$trim) {
    osl <- sapply(subsets, length)
    subsets <- lapply(subsets, function(sst) {
      if (!is.character(sst)) 
        colnames(alternative)[sst]
      else
        intersect(sst, colnames(alternative))
    })
  }
                                  
  # make sure that weights is a named list, compatible with colnames(alternative)
  if (many.weights) {
    names.weights <- names(weights)
    weights <- lapply(1:length(weights), function (i) {
      wt <- weights[[i]]
      if (!is.null(names(wt)))
        wt <- wt[names(wt) %in% colnames(alternative)]
      else 
        if (many.subsets && length(wt) == length(subsets[[i]]))
          names(wt) <- subsets[[i]]
        else if (length(wt) == ncol(alternative))
          names(wt) <- colnames(alternative)
      wt
    })
    names(weights) <- names.weights
    if (any(sapply(lapply(weights, names), is.null)))
      stop("weights input is not compatible with covariates input.")
  }
  
  # make subsets and weights compatible
  if (many.subsets && many.weights) {
    weights <- lapply(1:length(weights), function(i) {
      if (all(subsets[[i]] %in% names(weights[[i]])))
        weights[[i]][subsets[[i]]]
      else 
        stop("names of weights input incompatible with subsets input.")
    })
  }

  # trim zero weights
  if (many.weights) {
    if (any(unlist(weights)==0)) {
      weights <- lapply(weights, function(wt) wt[wt != 0])
      many.subsets <- FALSE   # to redo subsets
    }
  }
  
  # make subsets in case of short named weights
  if (many.weights && !many.subsets && any(sapply(weights, length) != ncol(alternative))) {
    subsets <- lapply(weights, names)
    many.subsets <- TRUE
  }
  
  # check missing values in subsets
  if (many.subsets && any(sapply(subsets, function(x) any(is.na(x))))) {
    stop("missing values in \"subsets\"")
  }
  
  # prepare progress info
  if (missing(trace)) trace <- gt.options()$trace && (many.weights || many.subsets)
  if (trace && (many.subsets || many.weights)) {
    if (many.subsets) 
      K <- length(subsets) 
    else 
      K <- length(weights) 
    digitsK <- trunc(log10(K))+1
  }

  # standardize and weight
  if (standardize || one.weight) {
    if (one.weight) {
      if (length(weights) != ncol(alternative)) 
        stop("length of \"weights\" does not match column count of \"alternative\"")
      all.weights <- weights
    } else
      all.weights <- rep(1, ncol(alternative))
    if (standardize) {
      vars <- apply(alternative,2,var) * (n-1)/n
      vars[vars == 0] <- 1
      all.weights <- all.weights / vars
    }
    alternative <- alternative * matrix(sqrt(all.weights), nrow(alternative), ncol(alternative), byrow=TRUE)
  }
                       
  # Make the test function
  funs <- switch(model,
    linear = .lineartest(response, Z=null, X=alternative, offset=offset, dir = directional, perms = permutations),
    logistic = .glmtest(response, Z=null, X=alternative, offset=offset, family = "binomial", dir = directional, perms=permutations),
    poisson = .glmtest(response, Z=null, X=alternative, offset=offset, family = "poisson", dir = directional, perms=permutations),
    multinomial = .multinomialtest(response, Z=null, X=alternative, dir = directional, perms=permutations),
    cox = .coxtest(response, Z=null, X=alternative, offset=offset, dir = directional, perms=permutations)
  )    
  test <- funs$test
                                                      
  # Do the test
  if ((!many.subsets) && (!many.weights)) {           # single weighting; single subset
    res <- test()
  } else {     
    L <- if (many.subsets) length(subsets) else length(weights)                 
    res <- sapply(1:L, function (i) {
      if (trace && L>1) {
        cat(rep("\b", 2*digitsK+3), i, " / ", K, sep="")
        flush.console()
      }
      if (!many.weights) {                                           # single weighting; many subsets
        uit <- test(subset=subsets[[i]]) 
      } else if (!many.subsets) {                                    # many weightings; single subset
        uit <- test(weights=weights[[i]]) 
      } else {                                                       # many weightings; many subsets
        uit <- test(subset=subsets[[i]], weights=weights[[i]])
      } 
      uit
    })
    if (many.subsets && !is.null(names(subsets)))
      colnames(res) <- names(subsets)
    else if (many.weights && !is.null(names(weights)))
      colnames(res) <- names(weights)
  }
  if (trace && (many.subsets || many.weights) && L>1) cat("\n")
  res <- t(res)
                                                    
  # return
  out <- new("gt.object")
  out@call <- call
  colnames(res) <- c("p-value", "Statistic", "Expected", "Std.dev", "#Cov")
  out@result <- res
  out@functions <- funs
  if (many.subsets) out@subsets <- subsets
  if (many.weights) out@weights <- weights
  out@directional <- directional > 0
  out@legend <- legend
  out@model <- model
  if (!is.null(alias)) alias(out) <- alias

  if (x) {
    out@null <- null
    out@alternative <- alternative
  }
  
  return(out)
}
