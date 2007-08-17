setClass("mlogit", representation(
    coefficients = "matrix",
    standard.err = "matrix",
    fitted.values = "matrix",
    x = "matrix",
    y = "matrix",
    formula = "formula",
    call = "call",
    df.null = "numeric",
    df.residual = "numeric",
    null.deviance = "numeric",
    deviance = "numeric",
    iter = "numeric",
    converged = "logical"
  )
)

setMethod("show", "mlogit", function(object) {
  cat("\n")
  cat(paste("Call:", deparse(object@call), "\n"))
  cat("\n")
  cat("Coefficients:\n")
  print(t(object@coefficients))
  cat("\n")
  cat(paste("Degrees of Freedom:", object@df.null, "Total (i.e. Null); ", object@df.residual, "Residual\n"))
  cat(paste("Null deviance:    ", round(object@null.deviance,2), "\n"))
  cat(paste("Residual deviance:", round(object@deviance,2), "\t"))
  g <- ncol(object@y)
  nullnulldf <- object@df.null + (g-1)
  AIC <- object@deviance + 2 * (nullnulldf - object@df.residual)
  cat(paste("AIC:", round(AIC,2) , "\n"))
})


setMethod("summary", "mlogit", function(object,...) {
  cat("\n")
  cat("Call:\n")
  cat(paste(deparse(object@call), "\n"))
  cat("\n")
  devres <- sqrt(-2 * (log(object@fitted.values)))
  devres[object@y==0] <- sqrt(-2 * (log(1-object@fitted.values)))[object@y==0]
  devres <- sign(object@y - object@fitted.values) * devres
  sumdevres <- t(round(apply(devres,2,quantile),4))
  colnames(sumdevres) <- c("Min", "1Q", "Median", "3Q", "Max")
  cat("Deviance Residuals:\n")
  show(sumdevres)
  cat("\n")
  cat("Coefficients:\n")
  g <- ncol(object@y)
  p <- ncol(object@x)
  for (ix in 1:g) {
    cat(paste("\nOutcome category ",colnames(object@coefficients)[ix], ":\n", sep = ""))
    coef <- object@coefficients[,ix]
    sterr <- object@standard.err[,ix]
    z <- coef / sterr
    p <- 2*pnorm(-abs(z))
    sumcoef <- cbind(coef, sterr, z, p)
    colnames(sumcoef) <- c("Estimate", "Std.Error", "z-value", "Pr(>|z|)")
    print(round(sumcoef, 4))
  }
  cat("\n")
  cat(paste("Null deviance:    ", round(object@null.deviance,3), "on", object@df.null, "degrees of freedom\n"))
  cat(paste("Residual deviance:", round(object@deviance,3), "on", object@df.residual, "degrees of freedom\n"))
  nullnulldf <- object@df.null + (g-1)
  AIC <- object@deviance + 2 * (nullnulldf - object@df.residual)
  cat(paste("AIC:", round(AIC,3) , "\n"))
  cat("\n")
  cat(paste("Number of Fisher Scoring iterations:", object@iter, "\n"))
})

setMethod("coefficients", "mlogit", function(object,...) {
  object@coefficients
})

setMethod("fitted.values", "mlogit", function(object,...) {
  object@fitted.values
})

setMethod("residuals", "mlogit", function(object,...) {
  object@y - object@fitted.values
})


mlogit <- function(formula, data, control = glm.control()) 
{
  # extract X from formula and data
  if (!is(formula, "formula")) stop("formula should be a formula object")
  if (!missing(data) && !is.data.frame(data)) stop("data should be a data.frame")

  # extract Y from formula and data  
  if (missing(data))
    data <- NULL
  Y <- factor(eval(formula[[2]], data))
  n <- length(Y)
  outs <- levels(Y)
  g <- length(outs)
  matrixY <- sapply(outs, function(out) { Y == out }) + 0
  
  dummyform <- as.formula(paste("rep(0,n) ~", formula[3]))
  dummyfit <- lm(dummyform, data = data, x = TRUE)
  X <- as.matrix(dummyfit$x)
  p <- ncol(X)
  
  # fit parameters
  # 1: starting values:
  gamma <- matrix(0,p,g)
  iter <- 0
  olddev <- -Inf
  finished <- FALSE
  
  # 2: Newton-Rhaphson algorithm
  while (!finished) {
    numerator <- sapply(1:g, function(out) exp(X %*% gamma[,out]))
    denominator <- rowSums(numerator)
    mu <- numerator / (denominator %o% rep(1,g))
    weightedX <- lapply(as.list(1:g), function(out) X * (mu[,out] %o% rep(1,p)))
    XWX <- matrix(0,p*g,p*g)
    for (ix in 1:g) {
      for (iy in 1:g) {
        if (ix == iy) {
          XWX[((ix-1)*p+1):(ix*p), ((ix-1)*p+1):(ix*p)] <- crossprod(weightedX[[ix]], X) - crossprod(weightedX[[ix]])
        } else {
          XWX[((ix-1)*p+1):(ix*p), ((iy-1)*p+1):(iy*p)] <- - crossprod(weightedX[[ix]], weightedX[[iy]])
        }
      }
    }
    eigs <- eigen(XWX, symmetric = TRUE)
    weightedEigens <- eigs$vectors[,1:(p*(g-1))] * (rep(1,p*g) %o% sqrt(1/eigs$values[1:(p*(g-1))]))
    XWX.MP <- crossprod(t(weightedEigens))
    gamma <- gamma + matrix(XWX.MP %*% as.vector(crossprod(X, matrixY - mu)), p, g)
    iter <- iter + 1
    dev <- -2 * sum(log(mu) * matrixY)  
    finished <- ( abs(dev - olddev) / (abs(dev) + 0.1) < control$epsilon ) | (iter >= control$maxit)
    olddev <- dev
  }
  if (iter == control$maxit) { warning("algorithm had not converged after ", control$maxit, " iterations") }    
  
  # 3: Calculate mu based on final estimates
  numerator <- sapply(1:g, function(out) exp(X %*% gamma[,out]))
  denominator <- rowSums(numerator)
  mu <- numerator / (denominator %o% rep(1,g))

  # Prepare output
  probs <- mu
  rownames(probs) <- names(Y)
  colnames(probs) <- outs
  coefs <- gamma
  rownames(coefs) <- colnames(X)
  colnames(coefs) <- outs
  call <- sys.call(0)
  nullmu <- rep(1,n) %o% colMeans(matrixY) 
  nulldev <- -2 * sum(log(nullmu) * matrixY)
  
  fit <- new("mlogit",
    coefficients = coefs,
    standard.err = matrix(sqrt(diag(XWX.MP)),p,g),
    x = X,
    y = matrixY,
    fitted.values = probs,
    call = call,
    formula = formula,
    df.null = (n-1) * (g-1),
    df.residual = (n-p) * (g-1),
    null.deviance = nulldev,
    deviance = dev,
    iter = iter,
    converged = (iter < control$maxit)
  )
  
  fit
}

.EQall <- function(X, IminH, fit) {
  mu <- fitted.values(fit)
  n <- nrow(mu)
  g <- ncol(mu)
  m <- nrow(X)
  range <- lapply(as.list(1:g), function(ix) ((ix-1)*n+1):(ix*n))
  matrixW <- matrix(0,g,n*g)
  for (ix in 1:g) 
    matrixW[,range[[ix]]] <- -t(mu)
  matrixW <- matrixW * (rep(1,g) %o% as.vector(mu))
  for (ix in 1:g) 
    matrixW[ix,range[[ix]]] <- matrixW[ix,range[[ix]]] + mu[,ix]
  XX <- crossprod(X) / m
  R <- matrix(0,n*g,n*g)
  for (ix in 1:g) {
    for(iy in 1:g) {
      R[range[[ix]],range[[iy]]] <- crossprod(IminH[range[[ix]],range[[iy]]], XX) %*% 
        IminH[range[[ix]],range[[iy]]]
    }
  }
  RW <- matrix(0,n*g,n*g)
  for (ix in 1:g) {
    for (iy in 1:g) {
      M <- matrix(0, n, n)
      for (iz in 1:g) {
        M <- M + R[range[[ix]], range[[iz]]] * (matrixW[iz, range[[iy]]] %o% rep(1,n))
      }
      RW[range[[ix]], range[[iy]]] <- M
    }
  }
  sum(diag(RW))
}
  
