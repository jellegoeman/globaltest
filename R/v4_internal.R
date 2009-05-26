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
# A small iterations utility
#==========================================================
iterations <- function(object) {
  if (is(object, "mlogit"))
    object@iter
  else if ("iter" %in% names(object))
    object$iter
  else 
    0
}    




#==========================================================
# Calculate the IminH matrix
#==========================================================
.getIminH <- function(fit) {
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
  } else if (is(fit, "mlogit")) {
    if (ncol(fit@x) == 1) {
      IminH <- NULL
    } else { 
      mu <- fitted.values(fit)
      p <- ncol(fit@x)
      g <- ncol(mu)
      n <- nrow(mu)
      weightedX <- lapply(as.list(1:g), function(out) fit@x * (mu[,out] %o% rep(1,p)))
      XWX <- matrix(0,p*g,p*g)
      for (ix in 1:g) {
        for (iy in 1:g) {
          if (ix == iy) {
            XWX[((ix-1)*p+1):(ix*p), ((ix-1)*p+1):(ix*p)] <- crossprod(weightedX[[ix]], fit@x) - crossprod(weightedX[[ix]])
          } else {
            XWX[((ix-1)*p+1):(ix*p), ((iy-1)*p+1):(iy*p)] <- - crossprod(weightedX[[ix]], weightedX[[iy]])
          }
        }
      }
      eigs <- eigen(XWX, symmetric = TRUE)
      weightedEigens <- eigs$vectors[,1:(p*(g-1))] * (rep(1,p*g) %o% sqrt(1/eigs$values[1:(p*(g-1))]))
      XWX.MP <- crossprod(t(weightedEigens))
      XXWXX <- matrix(0,n*g,n*g)
      for (ix in 1:g) {
        for (iy in 1:g) {
            XXWXX[((ix-1)*n+1):(ix*n), ((iy-1)*n+1):(iy*n)] <- 
              fit@x %*% XWX.MP[((ix-1)*p+1):(ix*p), ((iy-1)*p+1):(iy*p)] %*% t(fit@x)
        }
      }
      matrixW <- matrix(0,g,n*g)
      for (ix in 1:g) matrixW[,((ix-1)*n+1):(ix*n)] <- -t(mu)
      matrixW <- matrixW * (rep(1,g) %o% as.vector(mu))
      for (ix in 1:g) matrixW[ix,((ix-1)*n+1):(ix*n)] <- matrixW[ix,((ix-1)*n+1):(ix*n)] + mu[,ix]
      IminH <- matrix(0,n*g,n*g)
      for (ix in 1:g) {
        for (iy in 1:g) {
          M <- matrix(0, n, n)
          for (iz in 1:g) {
            M <- M + XXWXX[((iz-1)*n+1):(iz*n), ((iy-1)*n+1):(iy*n)] * (matrixW[ix, ((iz-1)*n+1):(iz*n)] %o% rep(1,n))
          }
          IminH[((ix-1)*n+1):(ix*n), ((iy-1)*n+1):(iy*n)] <- -M
        }
      }
      diag(IminH) <- diag(IminH) + 1
    }
  }
  IminH
}

#==========================================================
# A quick globaltest function for already prepared data
# This fuction chooses the appropriate globaltest function based 
# on the chosen model (asymptotic or gamma). 
#==========================================================
.globaltest <- function(gt, genesets, accuracy = 50) {
  if (!missing(genesets)) {
    gt@genesets <- genesets
  }
  model <- .model(gt)
  adjusted <- .adjusted(gt)
  if (gt@method == 2) {
    if (model == "linear") {
      res <- .linearglobaltestgamma(gt)
    } else if (model == "logistic") {
      if (adjusted) {
        res <- .adjustedlogisticglobaltestgamma(gt)
      } else {
        res <- .unadjustedlogisticglobaltestgamma(gt)
      }
    } else if (model == "survival") {
      if (adjusted) {
        res <- .adjustedsurvivalglobaltest(gt)
      } else {
        res <- .unadjustedsurvivalglobaltest(gt)
      }
    } else if (model == "multinomial") {
      if (adjusted) {
        res <- .adjustedmultinomialglobaltestgamma(gt)
      } else {
        res <- .unadjustedmultinomialglobaltestgamma(gt)
      }
    }
  } else if (gt@method == 3) {  
    if (model == "logistic") {
      res <- .logisticglobaltest(gt, accuracy)
    } else if (model == "linear") {
      res <- .linearglobaltest(gt, accuracy)
    } else  if (model == "survival") {
      if (adjusted) {
        res <- .adjustedsurvivalglobaltest(gt)
      } else {
        res <- .unadjustedsurvivalglobaltest(gt)
      }
    } else if (model == "multinomial") {
      if (adjusted) {
        res <- .adjustedmultinomialglobaltest(gt, accuracy)
      } else {
        res <- .unadjustedmultinomialglobaltest(gt, accuracy)
      }
    }
  }
  res
}


#==========================================================
# The globaltest for a linear model; asymptotic distribution
#==========================================================
.linearglobaltest <- function (gt, accuracy = 50) {
  Y <- fit(gt)$y - fitted.values(fit(gt))
  mu2 <- sum(Y*Y)/df.residual(fit(gt))
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
        lams <- eigen(R, symmetric = TRUE, only.values = TRUE)$values
      } else {
        if (.adjusted(gt)) {
          otherR <- crossprod(t(X %*% .IminH(gt))) / m
        } else {
          otherR <- crossprod(t(X)) / m
        }
        Q <- crossprod(X %*% Y) / (m * mu2)
        lams <- eigen(otherR, symmetric = TRUE, only.values = TRUE)$values
      }
      EQ <- sum(lams)
      varQ <- 2*sum(lams * lams)
      seQ <- sqrt(varQ)
      p.value <- .pAsymptotic(Q, .weed(lams, accuracy))
      out <- c(Q, EQ, seQ, p.value)
    }
    out
  }))
  res
}    

     
#==========================================================
# The globaltest for the logistic model; asymptotic distribution
#==========================================================
.logisticglobaltest <- function (gt, accuracy = 50) {
  Y <- fit(gt)$y - fitted.values(fit(gt))
  n <- length(Y)
  mu2 <- fit(gt)$weights
  res <- t(sapply(.genesets(gt), function(set) {
    X <- gt@eX[set,,drop=FALSE]
    m <- length(set)
    if (m == 0) {
      out <- rep(NA, 4)
    } else {
      if (m > n) {
        XX <- crossprod(X) / m
        if (.adjusted(gt)) {
          Q <- (crossprod(Y, XX) %*% Y)
          R <- crossprod(.IminH(gt), XX) %*% .IminH(gt)
          R <- R * (sqrt(mu2) %o% sqrt(mu2))
        } else {
          Q <- (crossprod(Y, XX) %*% Y) / mu2[1]
          R <- XX
        }
        lams <- eigen(R, symmetric = TRUE, only.values = TRUE)$values
      } else {
        if (.adjusted(gt)) {
          Q <- crossprod(X %*% Y) / m
          XIminH <- (X %*% .IminH(gt)) * (rep(1,m) %o% sqrt(mu2))
          otherR <- crossprod(t(XIminH)) / m
        } else {
          Q <- crossprod(X %*% Y) / (m * mu2[1])
          otherR <- crossprod(t(X)) / m
        }
        lams <- eigen(otherR, symmetric = TRUE, only.values = TRUE)$values
      }
      EQ <- sum(lams)
      varQ <- 2*sum(lams * lams)
      seQ <- sqrt(varQ)
      p.value <- .pAsymptotic(Q, .weed(lams, accuracy))
      out <- c(Q, EQ, seQ, p.value)
    }
    out
  }))
  res
}    


#==========================================================
# The globaltest for the unadjusted multinomial model; asymptotic distribution
#==========================================================
.unadjustedmultinomialglobaltest <- function (gt, accuracy = 50) {
  Y <- fit(gt)@y - fitted.values(fit(gt))
  g <- ncol(Y)
  n <- nrow(Y)
  mu <- fitted.values(fit(gt))[1,]
  matrixW <- -mu %o% mu + diag(mu)
  lamsW <- eigen(matrixW, symmetric = TRUE, only.values = TRUE)$values
  res <- t(sapply(.genesets(gt), function(set) {
    X <- gt@eX[set,,drop=FALSE]
    m <- length(set)
    if (m == 0) {
      out <- rep(NA, 4)
    } else {
      if (m > n) {
        XX <- crossprod(X) / m
        Q <- sum(sapply(1:g, function(out) crossprod(Y[,out], XX) %*% Y[,out]))
        lamsR <- eigen(XX, symmetric = TRUE, only.values = TRUE)$values
      } else {
        Q <- sum(sapply(1:g, function(out) crossprod(X %*% Y[,out]))) / m
        otherXX <- crossprod(t(X)) / m
        lamsR <- eigen(otherXX, symmetric = TRUE, only.values = TRUE)$values
      }
      lams <- lamsW %x% lamsR
      EQ <- sum(lams)
      varQ <- 2*sum(lams * lams)
      seQ <- sqrt(varQ)
      p.value <- .pAsymptotic(Q, .weed(lams, accuracy))
      out <- c(Q, EQ, seQ, p.value)
    }
    out
  }))
  res
}    


#==========================================================
# The globaltest for the adjusted multinomial model; asymptotic distribution
#==========================================================
.adjustedmultinomialglobaltest <- function (gt, accuracy = 50) {
  mu <- fitted.values(fit(gt))
  n <- nrow(mu)
  g <- ncol(mu)
  Y <- fit(gt)@y - mu
  range <- lapply(as.list(1:g), function(ix) ((ix-1)*n+1):(ix*n))
  matrixW <- matrix(0,g,n*g)
  for (ix in 1:g) 
    matrixW[,range[[ix]]] <- -t(mu)
  matrixW <- matrixW * (rep(1,g) %o% as.vector(mu))
  for (ix in 1:g) 
    matrixW[ix,range[[ix]]] <- matrixW[ix,range[[ix]]] + mu[,ix]
  res <- t(sapply(.genesets(gt), function(set) {
    X <- gt@eX[set,,drop=FALSE]
    m <- length(set)
    if (m == 0) {
      out <- rep(NA, 4)
    } else {
      IminH <- .IminH(gt)
      if (m > n) {
        XX <- crossprod(X) / m
        Q <- sum(sapply(1:g, function(out) crossprod(Y[,out], XX) %*% Y[,out]))
        R <- matrix(0,n*g,n*g)
        for (ix in 1:g) {
          for(iy in 1:g) {
            M <- matrix(0,n,n)
            for (iz in 1:g) {
              M <- M + crossprod(IminH[range[[iz]],range[[ix]]], XX) %*% 
                IminH[range[[iz]],range[[iy]]]
            }
            R[range[[ix]],range[[iy]]] <- M
          }
        }
        RW <- matrix(0,n*g,n*g)
        for (ix in 1:g) {
          for (iy in 1:g) {
            M <- matrix(0, n, n)
            for (iz in 1:g) {
              M <- M + R[range[[iz]], range[[iy]]] * (matrixW[ix, range[[iz]]] %o% rep(1,n))
            }
            RW[range[[ix]], range[[iy]]] <- M
          }
        }
        lams <- svd(RW, nu=0, nv=0)$d
      } else {
        Q <- sum(sapply(1:g, function(out) crossprod(X %*% Y[,out]))) / m
        XIminH <- matrix(0,m*g,n*g)
        for (ix in 1:g) {
          for (iy in 1:g) {
            XIminH[((ix-1)*m+1):(ix*m),range[[iy]]] <- X %*% IminH[range[[ix]], range[[iy]]]
          }
        }
        XIminHW <- matrix(0,m*g,n*g)
        for (ix in 1:g) {
          for (iy in 1:g) {
            M <- matrix(0, m, n)
            for (iz in 1:g) {
              M <- M + XIminH[((ix-1)*m+1):(ix*m), range[[iz]]] * (rep(1,m) %o% matrixW[iz, range[[iy]]])
            }
            XIminHW[((ix-1)*m+1):(ix*m), range[[iy]]] <- M
          }
        }        
        otherRW <- XIminHW %*% t(XIminH) / m
        lams <- eigen(otherRW, symmetric = TRUE, only.values = TRUE)$values
      }
      EQ <- sum(lams)
      varQ <- 2*sum(lams * lams)
      seQ <- sqrt(varQ)
      p.value <- .pAsymptotic(Q, .weed(lams, accuracy))
      out <- c(Q, EQ, seQ, p.value)
    }
    out
  }))
  res
}    




#==========================================================
# The globaltest for a linear model; gamma approximation
#==========================================================
.linearglobaltestgamma <- function (gt) {
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
# The globaltest for an unadjusted logistic model; gamma approximation
#==========================================================
.unadjustedlogisticglobaltestgamma <- function (gt) {
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
# The globaltest for an adjusted logistic model; gamma approximation
#==========================================================
.adjustedlogisticglobaltestgamma <- function (gt) {
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
# The globaltest for a multinomial model; gamma approximation
#==========================================================
.adjustedmultinomialglobaltestgamma <- function (gt) {
  mu <- fitted.values(fit(gt))
  n <- nrow(mu)
  g <- ncol(mu)
  kappai <- function(stuv) {
    case <- rowSums(outer(1:g, stuv, "=="))
    switch(sort(case, decreasing = TRUE)[1],
      { # sort(case) = c(1,1,1,1)
        muis <- mu[,case == 1][,1]
        muit <- mu[,case == 1][,2]
        muiu <- mu[,case == 1][,3]
        muiv <- mu[,case == 1][,4]
        -6*muis*muit*muiu*muiv
      },
      {
      switch(sort(case, decreasing = TRUE)[2],
        { # sort(case) = c(2,1,1,0) 
          muis <- mu[,case == 2]
          muit <- mu[,case == 1][,1]
          muiu <- mu[,case == 1][,2]
          2*muis*muit*muiu - 6*muis*muis*muit*muiu
        },
        { # sort(case) = c(2,2,0,0) 
          muis <- mu[,case == 2][,1]
          muit <- mu[,case == 2][,2]
          -muis*muit + 2*muis*muit*muit + 2*muis*muis*muit - 6*muis*muis*muit*muit
        })
      },
      { # sort(case) = c(3,1,0,0)
        muis <- mu[,case == 3]
        muit <- mu[,case == 1]
        -muis*muit + 6*muit*muis*muis - 6*muit*muis*muis*muis
      },
      { # sort(case) = c(4,0,0,0)
        muis <- mu[,case == 4]
        muis - 7*muis*muis + 12*muis*muis*muis - 6*muis*muis*muis*muis
      }
    )
  }
  Y <- fit(gt)@y - mu
  range <- lapply(as.list(1:g), function(ix) ((ix-1)*n+1):(ix*n))
  matrixW <- matrix(0,g,n*g)
  for (ix in 1:g) 
    matrixW[,range[[ix]]] <- -t(mu)
  matrixW <- matrixW * (rep(1,g) %o% as.vector(mu))
  for (ix in 1:g) 
    matrixW[ix,range[[ix]]] <- matrixW[ix,range[[ix]]] + mu[,ix]
  res <- t(sapply(.genesets(gt), function(tg) {
    X <- gt@eX[tg,,drop=FALSE]
    m <- length(tg)
    if (m == 0) {
      out <- rep(NA, 4)
    } else {
      IminH <- .IminH(gt)
      if (m > n) {
        XX <- crossprod(X) / m
        R <- matrix(0,n*g,n*g)
        for (ix in 1:g) {
          for(iy in 1:g) {
            M <- matrix(0,n,n)
            for (iz in 1:g) {
              M <- M + crossprod(IminH[range[[iz]],range[[ix]]], XX) %*% 
                IminH[range[[iz]],range[[iy]]]
            }
            R[range[[ix]],range[[iy]]] <- M
          }
        }
        Q <- sum(sapply(1:g, function(out) crossprod(Y[,out], XX) %*% Y[,out]))
        RW <- matrix(0,n*g,n*g)
        for (ix in 1:g) {
          for (iy in 1:g) {
            M <- matrix(0, n, n)
            for (iz in 1:g) {
              M <- M + R[range[[iz]], range[[iy]]] * (matrixW[ix, range[[iz]]] %o% rep(1,n))
            }
            RW[range[[ix]], range[[iy]]] <- M
          }
        }
        dRW <- sum(diag(RW))
        dRWRW <- sum(RW*RW)
        dRs <- lapply(as.list(1:g), function(out1) {
          lapply(as.list(1:g), function(out2) {
            diag(R[range[[out1]],range[[out2]]])
          })
        })
      } else { 
        Q <- sum(sapply(1:g, function(out) crossprod(X %*% Y[,out]))) / m
        XIminH <- matrix(0,m*g,n*g)
        for (ix in 1:g) {
          for (iy in 1:g) {
            XIminH[((ix-1)*m+1):(ix*m),range[[iy]]] <- X %*% IminH[range[[ix]], range[[iy]]]
          }
        }
        XIminHW <- matrix(0,m*g,n*g)
        for (ix in 1:g) {
          for (iy in 1:g) {
            M <- matrix(0, m, n)
            for (iz in 1:g) {
              M <- M + XIminH[((ix-1)*m+1):(ix*m), range[[iz]]] * (rep(1,m) %o% matrixW[iz, range[[iy]]])
            }
            XIminHW[((ix-1)*m+1):(ix*m), range[[iy]]] <- M
          }
        }        
        otherRW <- XIminHW %*% t(XIminH) / m
        dRW <- sum(diag(otherRW))
        dRWRW <- sum(otherRW*otherRW)
        XIminHXIminH <- XIminH*XIminH / m
        dRs <- lapply(as.list(1:g), function(ix) {
          lapply(as.list(1:g), function(iy) {
            colSums(XIminHXIminH[((ix-1)*m+1):(ix*m), range[[iy]],drop=FALSE])
          })
        })
      }
      EQ <- dRW
      varQ <- 2 * dRWRW + sum(sapply(1:g, function(out1) {
        sum(sapply(1:g, function(out2) {
          sum(sapply(1:g, function(out3) {
            sum(sapply(1:g, function(out4) {
              sum(dRs[[out1]][[out2]] * dRs[[out3]][[out4]] * kappai(c(out1, out2, out3, out4)))
            }))
          }))
        }))
      }))
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
# The globaltest for an unadjusted multinomial model; gamma approximation
#==========================================================
.unadjustedmultinomialglobaltestgamma <- function(gt) {
  Y <- fit(gt)@y - fitted.values(fit(gt))
  mu <- fitted.values(fit(gt))[1,]
  n <- nrow(Y)
  g <- ncol(Y)
  kappa <- function(stuv) {
    case <- rowSums(outer(1:g, stuv, "=="))
    switch(sort(case, decreasing = TRUE)[1],
      { # sort(case) = c(1,1,1,1)
        muis <- mu[case == 1][1]
        muit <- mu[case == 1][2]
        muiu <- mu[case == 1][3]
        muiv <- mu[case == 1][4]
        -6*muis*muit*muiu*muiv
      },
      {
      switch(sort(case, decreasing = TRUE)[2],
        { # sort(case) = c(2,1,1,0) 
          muis <- mu[case == 2]
          muit <- mu[case == 1][1]
          muiu <- mu[case == 1][2]
          2*muis*muit*muiu - 6*muis*muis*muit*muiu
        },
        { # sort(case) = c(2,2,0,0) 
          muis <- mu[case == 2][1]
          muit <- mu[case == 2][2]
          -muis*muit + 2*muis*muit*muit + 2*muis*muis*muit - 6*muis*muis*muit*muit
        })
      },
      { # sort(case) = c(3,1,0,0)
        muis <- mu[case == 3]
        muit <- mu[case == 1]
        -muis*muit + 6*muit*muis*muis - 6*muit*muis*muis*muis
      },
      { # sort(case) = c(4,0,0,0)
        muis <- mu[case == 4]
        muis - 7*muis*muis + 12*muis*muis*muis - 6*muis*muis*muis*muis
      }
    )
  }
  matrixW <- -mu %o% mu + diag(mu)
  res <- t(sapply(.genesets(gt), function(tg) {
    X <- gt@eX[tg,,drop=FALSE]
    m <- length(tg)
    if (m == 0) {
      out <- rep(NA, 4)
    } else {
      if (m > n) {
        XX <- crossprod(X) / m
        Q <- sum(sapply(1:g, function(out) crossprod(Y[,out], XX) %*% Y[,out]))
        traceXX <- sum(diag(XX))
        RW <- matrixW %x% XX
        dRW <- sum(diag(matrixW)) * sum(diag(XX))
        dRWRW <- sum(matrixW*matrixW) * sum(XX*XX)
      } else { 
        Q <- sum(sapply(1:g, function(out) crossprod(X %*% Y[,out]))) / m
        otherXX <- crossprod(t(X)) / m
        traceXX <- sum(diag(otherXX))
        dRW <- sum(diag(matrixW)) * sum(diag(otherXX))
        dRWRW <- sum(matrixW*matrixW) * sum(otherXX*otherXX)
      }
      EQ <- dRW
      varQ <- 2 * dRWRW + traceXX * sum(sapply(1:g, function(out1) {
        sum(sapply(1:g, function(out2) {
          sum(sapply(1:g, function(out3) {
            sum(sapply(1:g, function(out4) {
              matrixW[out1,out2] * matrixW[out3,out4] * kappa(c(out1, out2, out3, out4))
            }))
          }))
        }))
      }))
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
  

#==========================================================
# Removes extremely small eigenvalues 
#==========================================================
.weed <- function(lams, accuracy) {
  if (missing(accuracy)) {
    thresh <- 1/50
  } else {
    thresh <- 1/accuracy
  }
  lams <- -sort(-lams)
  m <- length(lams)
  while ((lams[1] > 0) && (lams[m] / lams[1] < thresh)) {
    q <- m-1
    r <- m-2
    lams[q] <- lams[q] + lams[m]
    while ((r > 0) && (lams[r] < lams[q])) {
      lams[r:q] <- mean(lams[r:q]) 
      r <- r - 1
    }
    m <- q
  }
  lams[1:m]
}

#==========================================================
# Calculates the asymptotic p-value using methods of
# Kotz, Johnson and Boyd (1967)
# Box (1954)
#==========================================================
.pAsymptotic <- function(x, lams, bet) {
  m <- length(lams)
  accuracy <- 10^-12
  if (lams[1] == 0) {
    upper <- 1
  } else {
    if (m == 1) {
        upper <- pchisq(x / lams, df = 1, lower.tail = FALSE)
    } else {
      # get the tuning parameter beta
      if (missing(bet)) {
        lams <- sort(lams)
        ruben <- 2 * lams[1] * lams[m] / (lams[1] + lams[m])
        harmonic <- 1/mean(1/lams)
        bet <- min(ruben, harmonic) * (1 - 10^-15)
      }
      # get an upper bound to the number of iterations needed
      A <- qnorm(.Machine$double.neg.eps)^2
      B <- x/bet
      maxiter <- trunc(0.5 * (A + B + sqrt(A*A + 2*A*B) - m))
      # starting values
      d <- numeric(maxiter)
      c <- numeric(maxiter+1)
      c[1] <- prod(sqrt(bet / lams))
      sumc <- c[1]
      chi <- pchisq(x / bet, df = m, lower.tail = FALSE) 
      partialsum <- c[1] * chi
      dbase <- (1 - bet /lams)
      ready <- FALSE
      mixture <- TRUE
      ix <- 1
      # iterate!
      while (!ready) {
        d[ix] <- 0.5 * sum(dbase^ix)
        c[ix+1] <- mean(c[1:ix] * d[ix:1])
        if (c[ix+1] < 0)
          mixture <- FALSE
        sumc <- sumc + c[ix+1]
        partialsum <- partialsum + c[ix+1] * chi
        chi <- pchisq(x / bet, df = m + 2 * ix + 2, lower.tail = FALSE)
        lower <- partialsum + (1 - sumc) * chi
        upper <- partialsum + 1 - sumc
        if (mixture)
          ready <- ((upper - lower) / (upper + lower) < 10^-5) || (ix == maxiter) || (upper < accuracy)
        else {
          ready <- TRUE
          upper <- .pAsymptotic(x, lams, mean(c(bet, min(lams))))
        }
        ix <- ix + 1
      }
    }
  }
  if (upper < accuracy)
    upper <- 0
  upper
}
