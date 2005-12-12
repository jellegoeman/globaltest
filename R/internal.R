#==========================================================
# .onLoad is called when the package is loaded
#==========================================================
.onLoad <- function(lib, pkg) require(methods)

.onAttach <- function(lib, pkg) {
    if(interactive() && .Platform$OS.type == "windows" && .Platform$GUI == "Rgui"){
      addVigs2WinMenu("globaltest")
    }
}
#==========================================================


#==========================================================
# Calculate the IminH matrix
#==========================================================
getIminH <- function(fit) {
  if (is(fit, "glm")) {
    if (fit$rank == 1) {
      IminH <- NULL
    } else {
      matrixW <- diag(fit$weights)
      IminH <- - matrixW %*% (fit$x %*% solve(crossprod(fit$x, matrixW) %*% fit$x, t(fit$x)))
      diag(IminH) <- diag(IminH) + 1
    }
  } else if (is(fit, "lm")) {
    if (fit$rank == 1) {
      IminH <- NULL
    } else { 
      IminH <- - fit$x %*% solve(crossprod(fit$x), t(fit$x))
      diag(IminH) <- diag(IminH) + 1
    }
  } else if (is(fit, "coxph.null")) {
    IminH <- NULL
  } else if (is(fit, "coxph")) {
    expci <- exp(fit$linear.predictors)
    n <- length(expci)
    times <- as.vector(fit$y)[1:n]
    d <- as.vector(fit$y)[(n+1):(2*n)] 
    dtimes <- unique(times[d == 1]) 
    nd <- length(dtimes)
    matrixO <- outer(times, dtimes, "==") * matrix(d, n, nd)
    atrisk <- outer(times, dtimes, ">=")
    hazinc <- as.vector(1 / (expci %*% atrisk))
    matrixP <- outer(expci, hazinc, "*") * atrisk
    matrixPO <- matrixP %*% t(matrixO)
    matrixW <- diag(rowSums(matrixPO)) - matrixPO %*% t(matrixPO)
    IminH <- - matrixW %*% (fit$x %*% solve(crossprod(fit$x, matrixW) %*% fit$x, t(fit$x)))
    diag(IminH) <- diag(IminH) + 1
  }
  IminH
}



#==========================================================
# The globaltest for a linear model
#==========================================================
.linearglobaltest <- function (gt) {
  Y <- fit(gt)$y - fitted.values(fit(gt))
  mu2 <- sum(Y*Y)/df.residual(fit(gt))
  nadjust <- ncol(fit(gt)$x)
  n <- length(Y)
  res <- t(sapply(.genesets(gt), function(tg) {
    X <- gt@eX[tg,,drop=FALSE]
    m <- length(tg)
    if (m == 0) {
      out <- rep(NA, 4)
    } else {
      if (m > n) {
        XX <- crossprod(X) / m
        if (.adjusted(gt)) {
          R <- crossprod(.IminH(gt), XX) %*% .IminH(gt)
        } else {
          R <- XX
        }
        Q <- (crossprod(Y, XX) %*% Y) / mu2
        trR <- sum(diag(R))
        trRR <- sum(R*R)          # trace(R^2) 
        tr2R <- trR * trR         # trace(R)^2
      } else {
        if (.adjusted(gt)) {
          otherR <- crossprod(t(X %*% .IminH(gt))) / m
        } else {
          otherR <- crossprod(t(X)) / m
        }
        Q <- crossprod(X %*% Y) / (m * mu2)
        trR <- sum(diag(otherR))
        trRR <- sum(otherR * otherR)    # trace(R^2) 
        tr2R <- trR * trR               # trace(R)^2
      }
      EQ <- trR
      varQ <- (2 / (n - nadjust + 2)) * ( (n - nadjust) * trRR - tr2R )
      seQ <- sqrt(varQ)
      scl <- varQ / (2 * EQ)
      dfr <- EQ / scl
      p.value <- pchisq(Q / scl, dfr, lower.tail = FALSE)
      out <- c(Q, EQ, seQ, p.value)
    }
    out
  }))
  res
}    

#==========================================================
# The globaltest for an unadjusted logistic model
#==========================================================
.unadjustedlogisticglobaltest <- function (gt) {
  Y <- fit(gt)$y - fitted.values(fit(gt))
  n <- length(Y)
  mu <- fitted.values(fit(gt))[1]
  mu2 <- fit(gt)$weights[1]
  K <- (1 - 6 * mu + 6 * mu^2 ) / mu2
  res <- t(sapply(.genesets(gt), function(set) {
    X <- gt@eX[set,,drop=FALSE]
    m <- length(set)
    if (m == 0) {
      out <- rep(NA, 4)
    } else {
      if (m > n) {
        XX <- crossprod(X) / m
        Q <- (crossprod(Y, XX) %*% Y) / mu2
        diagR <- colSums(X * X) / m
        trR <- sum(diagR)
        trRR <- sum(XX * XX)        # trace(R^2) 
        tr2R <- trR * trR           # trace(R)^2
        trR2 <- sum(diagR * diagR)  # trace(R.^2)
      } else {
        otherR <- crossprod(t(X)) / m
        Q <- crossprod(X %*% Y) / (m * mu2)
        trR <- sum(diag(otherR))
        trRR <- sum(otherR * otherR)    # trace(R^2) 
        tr2R <- trR * trR               # trace(R)^2
        diagR <- colSums(X * X) / m
        trR2 <- sum(diagR * diagR)      # trace(R.^2)
      }
      EQ <- trR
      varQ <- K * ( trR2 - tr2R / n ) + 2 * trRR - 2 * tr2R / (n-1)
      seQ <- sqrt(varQ)
      scl <- varQ / (2 * EQ)
      dfr <- EQ / scl
      p.value <- pchisq(Q / scl, dfr, lower.tail = FALSE)
      out <- c(Q, EQ, seQ, p.value)
    }
    out
  }))
  res
}    


#==========================================================
# The globaltest for an adjusted logistic model
#==========================================================
.adjustedlogisticglobaltest <- function (gt) {
  Y <- fit(gt)$y - fitted.values(fit(gt))
  n <- length(Y)
  mu <- fitted.values(fit(gt))
  mu2 <- fit(gt)$weights
  mu4 <- mu * (1-mu)^4 + (1-mu) * mu^4
  kurt <- mu4 - 3 * mu2 * mu2
  res <- t(sapply(.genesets(gt), function(set) {
    X <- gt@eX[set,,drop=FALSE]
    m <- length(set)
    if (m == 0) {
      out <- rep(NA, 4)
    } else {
      if (m > n) {
        XX <- crossprod(X) / m
        Q <- (crossprod(Y, XX) %*% Y)
        R <- crossprod(.IminH(gt), XX) %*% .IminH(gt)
        RV <- R * (rep(1,n) %o% mu2)
        diagR <- diag(R)
        traceRV <- sum(diag(RV))
        traceRVRV <- sum(RV * t(RV))
      } else {
        Q <- crossprod(X %*% Y) / m
        XIminH <- X %*% .IminH(gt)
        diagR <- colSums(XIminH * XIminH) / m
        XVX <- crossprod(t(XIminH * (rep(1,m) %o% sqrt(mu2)))) / m
        traceRV <- sum(diag(XVX))
        traceRVRV <- sum(XVX * XVX)
      }
      EQ <- traceRV
      varQ <-  sum(diagR * diagR * kurt) + 2 * traceRVRV
      seQ <- sqrt(varQ)
      scl <- varQ / (2 * EQ)
      dfr <- EQ / scl
      p.value <- pchisq(Q / scl, dfr, lower.tail = FALSE)
      out <- c(Q, EQ, seQ, p.value)
    }
    out
  }))
  res
}    
       
#==========================================================
# The  globaltest for an unadjusted survival model
# This function is far from optimal in terms of speed
# TODO: optimize!
#==========================================================
.unadjustedsurvivalglobaltest <- function (gt) {
  expci <- exp(fit(gt)$linear.predictors)
  n <- length(expci)
  times <- as.vector(fit(gt)$y)[1:n]
  d <- as.vector(fit(gt)$y)[(n+1):(2*n)] 
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
  Y <- rowSums(matrixPO) - d
  res <- t(sapply(.genesets(gt), function(tg) {
    X <- gt@eX[tg,,drop=FALSE]
    m <- length(tg)
    if (m == 0) {
      out <- rep(NA, 4)
    } else {
      XX <- crossprod(X) / m
      R <- XX
      Q <- crossprod(Y, XX) %*% Y
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
      EQ <- sum(diag(R %*% matrixW))   
      tussen <- diag(R) %o% rep(1,ncol(matrixP)) + 2 * R %*% (matrixM - matrixP)
      matrixT <- tussen - rep(1,n) %o% colSums(matrixP * tussen)    
      varQ <- sum(matrixPO * ((matrixT * matrixT) %*% t(matrixO)))
      Q <- (Q - EQ) / sqrt(varQ)
      p.value <- pnorm(-Q)
      EQ <- 0
      seQ <- 1
      out <- c(Q, EQ, seQ, p.value)
    }
    out
  }))
  res
}

#==========================================================
# The  globaltest for an unadjusted survival model
# This function is optimized for m >= n
# TODO: optimize for m < n
#==========================================================
.adjustedsurvivalglobaltest <- function (gt) {
  expci <- exp(fit(gt)$linear.predictors)
  n <- length(expci)
  times <- as.vector(fit(gt)$y)[1:n]
  d <- as.vector(fit(gt)$y)[(n+1):(2*n)] 
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
  Y <- rowSums(matrixPO) - d
  res <- t(sapply(.genesets(gt), function(tg) {
    X <- gt@eX[tg,,drop=FALSE]
    m <- length(tg)
    if (m == 0) {
      out <- rep(NA, 4)
    } else {
      XX <- crossprod(X) / m
      R <- crossprod(.IminH(gt), XX) %*% .IminH(gt)
      Q <- crossprod(Y, XX) %*% Y
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
      EQ <- sum(diag(R %*% matrixW))   
      tussen <- diag(R) %o% rep(1,ncol(matrixP)) + 2 * R %*% (matrixM - matrixP)
      matrixT <- tussen - rep(1,n) %o% colSums(matrixP * tussen)    
      varQ <- sum(matrixPO * ((matrixT * matrixT) %*% t(matrixO)))
      Q <- (Q - EQ) / sqrt(varQ)
      p.value <- pnorm(-Q)
      EQ <- 0
      seQ <- 1
      out <- c(Q, EQ, seQ, p.value)
    }
    out
  }))
  res
}
    
