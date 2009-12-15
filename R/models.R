##############################################
# Linear model (normal errors)
##############################################
.lineartest <- function(response, Z, X, offset, dir, perms) {

  # establish dimensions of the problem
  m <- ncol(Z)
  n <- nrow(Z)

  # fit the null model
  form <- formula(.makeFormula(m>0, !is.null(offset)))
  Y <- residuals(lm(form))
  sumYY <- sum(Y*Y)

  # adjust the alternative
  if (m > 0)
    X <- X - Z %*% solve(crossprod(Z), crossprod(Z, X))
  # set columns numerically zero to exactly zero
  csm <- colSums(X*X)
  X[,csm < max(csm)*1e-14] <- 0

  if (perms > 0) {
    # all permutations if possible
    if (npermutations(Y) <= perms) {
      random <- FALSE
      permY <- allpermutations(Y)
    } else {
      # otherwise random permutations
      random <- TRUE
      permY <- replicate(perms, Y[sample(n)])
    }
  }

  # Calculate the test statistic
  # 1: asymptotic version
  if (perms == 0) {
    test <- function(subset, weights, calculateP = TRUE) {
      if (!missing(subset))
        X <- X[,subset,drop=FALSE]
      p <- ncol(X)
      if (p == 0) {
        return(c(p = NA, S = NA, ES = NA, sdS = NA, ncov=0))
      } else {
        if(!missing(weights))
          X <- X * matrix(ifelse(weights>0, sqrt(weights), -sqrt(weights)), n, p, byrow = TRUE)
  
        if (p > n) {
          XX <- crossprod(t(X))
          if (dir) XX <- XX + dir * outer(rowSums(X), rowSums(X))
          norm.const <- sum(diag(XX)) / 100
          S <- sum(Y * (XX %*% Y)) / sum(Y*Y)
          lams <- eigen(XX, symmetric = TRUE, only.values=TRUE)$values
          lams[1:(n-m)] <- lams[1:(n-m)] - S
  
          # term needed for mean and variance of S for z-score
          var.num <- 2*sum(XX*XX)
  
        } else {
          if (dir && (p > 1)) X <- cbind(X, sqrt(dir) * rowSums(X))
  
          norm.const <- sum(X * X) / 100
          xy <- crossprod(X, Y)
          S <- sum(xy * xy) / sumYY
          lams <- eigen(crossprod(X), symmetric = TRUE, only.values=TRUE)$values
          if (length(lams) < n) lams <- c(lams, numeric(n-length(lams)))
          lams[1:(n-m)] <- lams[1:(n-m)] - S
  
          # term needed for mean and variance of S for z-score
          tr.term <- crossprod(X)
          var.num <- 2*sum(tr.term*tr.term)
        }
  
        # mean and variance of S for z-score
        # use series approximation as given in Paolella (2003) 319:
        mu.num <- sum(lams) + (n-m)*S
        mu.den <- n-m
        var.den <- 2*(n-m)
        cov.term <- 2*mu.num
        ES <- (mu.num/mu.den) * (1 - cov.term/(mu.num*mu.den) + var.den/(mu.den^2))
        VarS <- (mu.num^2/mu.den^2) * (var.num/(mu.num^2)  + var.den/(mu.den^2) - 2*cov.term/(mu.num*mu.den))
        
        # repair for case all(X == 0)
        if (norm.const == 0) {
          norm.const <- 1
          ES <- 0
          VarS <- 0
        }
                          
        # calculate the p-value
        if (calculateP)
          p.value <- .getP(lams)
        else
          p.value <- NA
                                        
        # give back
        return(c(p = p.value, S = S/norm.const, ES = ES/norm.const, sdS = sqrt(VarS)/norm.const, ncov=p))
      }
    }
  } else {
  # 2: permutation version
    test <- function(subset, weights, calculateP = TRUE) {          # calculateP not used in permutation testing
      if (!missing(subset))
        X <- X[,subset,drop=FALSE]
      p <- ncol(X)
      if (p == 0) {
        return(c(p = NA, S = NA, ES = NA, sdS = NA, ncov=0))
      } else {
        if(!missing(weights))
          X <- X * matrix(ifelse(weights>0, sqrt(weights), -sqrt(weights)), n, p, byrow = TRUE)
  
        if (p > n) {
  
          XX <- crossprod(t(X))
          if (dir) XX <- XX + dir * outer(rowSums(X), rowSums(X))
  
          norm.const <- sum(diag(XX)) * sumYY / 100
  
          S <- sum(Y * (XX %*% Y))
          permS <- colSums(permY * (XX %*% permY))
  
        } else {
          if (dir && (p > 1)) X <- cbind(X, sqrt(dir) * rowSums(X))
  
          norm.const <- sum(X * X) * sumYY / 100
          xy <- crossprod(X, Y)
          S <- sum(xy * xy)
  
          permxy <- crossprod(X, permY)
          permS <- colSums(permxy * permxy)
  
        }
  
        # mean and variance of S for z-score
        ES <- mean(permS)
        VarS <- var(permS)

        # repair for case all(X == 0)
        if (norm.const == 0) 
          norm.const <- 1
  
        # calculate the p-value
        p.value <- mean(S <= permS)
  
        # give back
        return(c(p = p.value, S = S/norm.const, ES = ES/norm.const, sdS = sqrt(VarS)/norm.const, ncov = p))
      }
    }
  }

  # Get the adjusted X matrix
  getX <- function(subset, weights) {
    if (!missing(subset))
      X <- X[,subset,drop=FALSE]
    p <- ncol(X)
    if(!missing(weights))
      X <- X * matrix(ifelse(weights>0, sqrt(weights), -sqrt(weights)), n, p, byrow = TRUE)
    X
  }

  # Subjectsplot function
  subjectplot <- function(subset, weights, mirror) {
    if (!missing(subset))
      X <- X[,subset,drop=FALSE]
    p <- ncol(X)
    if(!missing(weights))
      X <- X * matrix(ifelse(weights>0, sqrt(weights), -sqrt(weights)), n, p, byrow = TRUE)

    XX <- crossprod(t(X))
    if (dir) XX <- XX + dir * outer(rowSums(X), rowSums(X))
    if (mirror) {
      R <- XX %*% Y * sign(Y)
      ER <- diag(XX) * Y * sign(Y)
    } else {
      R <- XX %*% Y
      ER <- diag(XX) * Y
    }

    norm.const <- mean(abs(ER))

    XXmd <- XX
    diag(XXmd) <- 0
    if (m > 0) XXmd <- XXmd - (XXmd %*% Z) %*% solve(crossprod(Z), t(Z))
    varR <- colSums(XXmd * XXmd) * sumYY / (n-m)

    cbind(R=R/norm.const, ER=ER/norm.const, sdR=sqrt(varR)/norm.const, Y=Y)
  }

  # permutations for permutations histogram
  permutations <- function(subset, weights) {
      if (!perms)
        stop("histogram only for permutation version of the test")

      if (!missing(subset))
        X <- X[,subset,drop=FALSE]
      p <- ncol(X)
      if(!missing(weights))
        X <- X * matrix(ifelse(weights>0, sqrt(weights), -sqrt(weights)), n, p, byrow = TRUE)

      if (p > n) {

        XX <- crossprod(t(X))
        if (dir) XX <- XX + dir * outer(rowSums(X), rowSums(X))

        norm.const <- sum(diag(XX)) * sumYY / 100

        S <- sum(Y * (XX %*% Y))
        permS <- colSums(permY * (XX %*% permY))

      } else {
        if (dir && (p > 1)) X <- cbind(X, sqrt(dir) * rowSums(X))

        xy <- crossprod(X, Y)
        S <- sum(xy * xy)
        norm.const <- sum(X * X) * sumYY / 100

        permxy <- crossprod(X, permY)
        permS <- colSums(permxy * permxy)

      }

      return(list(S = S/norm.const, permS = permS/norm.const))
  }

  df <- function() {
    c(n,m,ncol(X))
  }

  nperms <- function() {
    if (perms)
      c(n = ncol(permY), random = random)
    else
      c(0, FALSE)
  }

  cov.names <- function(subset)
    if (missing(subset))
      colnames(X)
    else
      colnames(X[,subset,drop=FALSE])

  dist.cor <- function(subset, weights, transpose = FALSE) {
    if (!missing(subset))
      X <- X[,subset,drop=FALSE]
    p <- ncol(X)
    if(!missing(weights))
      X <- X * matrix(ifelse(weights>0, sqrt(weights), -sqrt(weights)), n, p, byrow = TRUE)
    if (transpose) 
      X <- t(X)
    sds <- sqrt(colSums(X*X))
    sds[sds < max(sds)*1e-14] <- Inf
    crossprod(X) / outer(sds,sds)
  }
  
  getY <- function() Y

  positive <- function(subset) {
    if (!missing(subset))
      X <- X[,subset,drop=FALSE] 
    xy <- crossprod(X, Y)
    drop(xy >= 0)+0
  }

  return(list(test=test,
    getX = getX,
    cov.names = cov.names,
    cor=dist.cor,
    getY = getY,
    subjectplot = subjectplot,
    positive = positive,
    permutations = permutations,
    nperms = nperms,
    df=df))

}


##############################################
# Generalized linear model (canonical link functions only)
##############################################
.glmtest <- function(response, Z, X, offset, family, dir, perms) {

  # establish dimensions of the problem
  m <- ncol(Z)
  n <- nrow(Z)

  # fit the null model
  form <- formula(.makeFormula(m>0, !is.null(offset)))
  null.fit <- glm(form, family = family)

  Y <- null.fit$y - fitted.values(null.fit)
  W <- null.fit$weights          
  if (all(W==W[1]))
    out <- .lineartest(as.numeric(response), Z, X, offset, dir, perms)
  else {
    if ((perms > 0) && (!all(W == W[1])))
      stop("Generalized linear model with covariates:\n\tPermutation testing is not possible.", call.=FALSE)
  
    if (perms == 0) {
      sqrtW <- sqrt(W)
      # adjust the alternative
      if (m > 0) {
        ZWhf <- Z * matrix(sqrtW, n, m)
        ZW <- Z * matrix(W, n, m)
        X <- X - Z %*% solve(crossprod(ZWhf), t(ZW)) %*% X
        csm <- colSums(X*X)
        X[,csm < max(csm)*1e-14] <- 0
        ZWZinvZW <- solve(crossprod(ZWhf), t(ZWhf))
      }
    } else {
      equal.size.groups <- all(abs(Y) == abs(Y[1])) && (sum(Y) == 0)
      if (equal.size.groups)
        nperms <- npermutations(Y[-1])
      else
        nperms <- npermutations(Y)
      # all permutations if possible
      if (nperms <= perms) {
        random <- FALSE
        if (equal.size.groups)
          permY <- rbind(Y[1], allpermutations(Y[-1]))
        else
          permY <- allpermutations(Y)
      } else {
        # otherwise random permutations
        random <- TRUE
        permY <- replicate(perms, Y[sample(n)])
      }
      # adjust the alternative
      if (m > 0)
        X <- X - matrix(colMeans(X), n, ncol(X), byrow=TRUE)
    }
  
    norm.const <- (n-m) / 100
  
    # calculate the test statistic
    # 1: asymptotic version
    if (perms == 0) {
      test <- function(subset, weights, calculateP = TRUE) {
        if (!missing(subset))
          X <- X[,subset,drop=FALSE]
        p <- ncol(X)
        if (p == 0) {
          return(c(p = NA, S = NA, ES = NA, sdS = NA, ncov=0))
        } else {
          if(!missing(weights))
            X <- X * matrix(ifelse(weights>0, sqrt(weights), -sqrt(weights)), n, p, byrow = TRUE)
    
          XX <- crossprod(t(X))
          if (dir) XX <- XX + dir * outer(rowSums(X), rowSums(X))
                         
          S <- sum(Y * (XX %*% Y)) / sum(Y*Y*diag(XX))
    
          # get lambdas
          XX <- XX * outer(sqrtW, sqrtW)
          dXX <- diag(diag(XX))
          if (m > 0) {
            dXX <- dXX - (dXX %*% ZWhf) %*% ZWZinvZW
            dXX <- dXX - crossprod(ZWZinvZW, crossprod(ZWhf, dXX))
          }
    
          if (is.nan(S))
            lams <- 0
          else
            lams <- eigen(XX - S*dXX, symmetric = TRUE, only.values=TRUE)$values
  
          # calculate the p-value
          if (calculateP)
            p.value <- .getP(lams)
          else
            p.value <- NA
    
          # mean and variance of S for z-score                                                           
          # use series approximation as given in Paolella (2003) 319:
          mu.num <- sum(diag(XX))
          mu.den <- sum(diag(dXX))
          var.num <- 2*sum(XX*XX)
          var.den <- 2*sum(dXX*dXX)
          cov.term <- 2*sum(XX*dXX)
    
          ES <- (mu.num/mu.den) * (1 - cov.term/(mu.num*mu.den) + var.den/(mu.den^2))
          VarS <- (mu.num^2/mu.den^2) * (var.num/(mu.num^2)  + var.den/(mu.den^2) - 2*cov.term/(mu.num*mu.den))
    
          # repair for case all(X == 0)
          if (is.na(S)) {
            S <- 0
            ES <- 0
            VarS <- 0
            p.value <- 1
          }
    
          # give back
          return(c(p = p.value, S = S/norm.const, ES = ES/norm.const, sdS = sqrt(VarS)/norm.const, ncov=p))
        }
      }
    } else {
    # 2: permutation version
      test <- function(subset, weights, calculateP=TRUE) {  # calculateP not used in permutation testing
        if (!missing(subset))
          X <- X[,subset,drop=FALSE]
        p <- ncol(X)
        if (p == 0) {
          return(c(p = NA, S = NA, ES = NA, sdS = NA, ncov=0))
        } else {
          if(!missing(weights))
            X <- X * matrix(ifelse(weights>0, sqrt(weights), -sqrt(weights)), n, p, byrow = TRUE)
    
          if (p > n) {
            XX <- crossprod(t(X))
            if (dir) XX <- XX + dir * outer(rowSums(X), rowSums(X))
    
            S <- sum(Y * (XX %*% Y)) / sum(Y*Y*diag(XX))
            permS <- colSums(permY * (XX %*% permY)) / drop(diag(XX) %*% (permY * permY))
    
          } else {
            if (dir) X <- X + dir * matrix(rowSums(X), nrow(X), ncol(X))
    
            diagXX <- rowSums(X * X)
            xy <- crossprod(X, Y)
            S <- sum(xy * xy) / sum(Y*Y*diagXX)
    
            permxy <- crossprod(X, permY)
            permS <- colSums(permxy * permxy) / drop(diagXX %*% (permY * permY))
    
          }
    
          # mean and variance of S for z-score
          ES <- mean(permS)
          VarS <- var(permS)
    
          # calculate the p-value
          p.value <- mean(S <= permS)
    
          # give back
          return(c(p = p.value, S = S/norm.const, ES = ES/norm.const, sdS = sqrt(VarS)/norm.const, ncov=p))
        }
      }
    }
  
    # Get the adjusted X matrix
    getX <- function(subset, weights) {
      if (!missing(subset))
        X <- X[,subset,drop=FALSE]
      p <- ncol(X)
      if(!missing(weights))
        X <- X * matrix(ifelse(weights>0, sqrt(weights), -sqrt(weights)), n, p, byrow = TRUE)
      X
    }
  
    # Subjectsplot function
    subjectplot <- function(subset, weights, mirror) {
      if (!missing(subset))
        X <- X[,subset,drop=FALSE]
      p <- ncol(X)
      if(!missing(weights))
        X <- X * matrix(ifelse(weights>0, sqrt(weights), -sqrt(weights)), n, p, byrow = TRUE)
  
      XX <- crossprod(t(X))
      if (dir) XX <- XX + dir * outer(rowSums(X), rowSums(X))
      if (mirror) {
        R <- XX %*% Y * sign(Y)
        ER <- diag(XX) * Y * sign(Y)
      } else {
        R <- XX %*% Y
        ER <- diag(XX) * Y
      }
  
      norm.const <- mean(abs(ER))
  
      XXmd <- XX
      diag(XXmd) <- 0
      if (m > 0) XXmd <- XXmd - (XXmd %*% ZW %*% solve(crossprod(ZWhf), t(Z)))
      varR <- drop((XXmd * XXmd) %*% W)
  
      cbind(R=R/norm.const, ER=ER/norm.const, sdR=sqrt(varR)/norm.const, Y=Y)
    }
  
    # permutations for permutations histogram
    permutations <- function(subset, weights) {
        if (!perms)
          stop("histogram only for permutation version of the test")
  
        if (!missing(subset))
          X <- X[,subset,drop=FALSE]
        p <- ncol(X)
        if(!missing(weights))
          X <- X * matrix(ifelse(weights>0, sqrt(weights), -sqrt(weights)), n, p, byrow = TRUE)
  
        if (p > n) {
          XX <- crossprod(t(X))
          if (dir) XX <- XX + dir * outer(rowSums(X), rowSums(X))
  
          S <- sum(Y * (XX %*% Y)) / sum(Y*Y*diag(XX))
          permS <- colSums(permY * (XX %*% permY)) / drop(diag(XX) %*% (permY * permY))
  
        } else {
          if (dir) X <- X + dir * matrix(rowSums(X), nrow(X), ncol(X))
  
          diagXX <- rowSums(X * X)
          xy <- crossprod(X, Y)
          S <- sum(xy * xy) / sum(Y*Y*diagXX)
  
          permxy <- crossprod(X, permY)
          permS <- colSums(permxy * permxy) / drop(diagXX %*% (permY * permY))
        }
  
        return(list(S = S/norm.const, permS = permS/norm.const))
    }
  
    df <- function() {
      c(n,m,ncol(X))
    }
  
    nperms <- function() {
      if (perms)
        c(n = ncol(permY), random = random)
      else
        c(0, FALSE)
    }
  
    cov.names <- function(subset)
      if (missing(subset))
        colnames(X)
      else
        colnames(X[,subset,drop=FALSE])
  
    dist.cor <- function(subset, weights, transpose = FALSE) {
      if (!missing(subset))
        X <- X[,subset,drop=FALSE]
      p <- ncol(X)
      if(!missing(weights))
        X <- X * matrix(ifelse(weights>0, sqrt(weights), -sqrt(weights)), n, p, byrow = TRUE)
      if (transpose) 
        X <- t(X)
      sds <- sqrt(colSums(X*X))
      sds[sds < max(sds)*1e-14] <- Inf
      crossprod(X) / outer(sds,sds)
    }
    
    getY <- function() Y
  
    positive <- function(subset, weights) {
      if (!missing(subset))
        X <- X[,subset,drop=FALSE] 
      xy <- crossprod(X, Y)
      drop(xy >= 0)+0
    }
  
    out <- list(test=test,
      getX = getX,
      cov.names = cov.names,
      cor=dist.cor,
      getY = getY,
      subjectplot = subjectplot,
      permutations = permutations,
      nperms = nperms,
      positive = positive,
      df=df)
  }
  
  return(out)
}


##############################################
# Cox proportional hazards model
##############################################
.coxtest <- function (response, Z, X, offset, dir, perms) {
                                 
  # fit the model
  m <- ncol(Z)
  form <- formula(.makeFormula(m>0, !is.null(offset)))
  fit <- coxph(form, method = "breslow")
  expci <- exp(fit$linear.predictors)
  n <- length(expci)

  # check possibility of permutations
  if ((perms > 0) && (m > 0))
    stop("Cox model with covariates:\n\tPermutation testing is not possible.", call.=FALSE)

  # extract survival times
  times <- as.vector(fit$y)[1:n]
  d <- as.vector(fit$y)[(n+1):(2*n)]
  dtimes <- unique(times[d == 1])
  nd <- length(dtimes)

  # any ties present?
  matrixO <- outer(times, dtimes, "==") * matrix(d, n, nd)
  ties <- any(colSums(matrixO) > 1)

  # construct relevant matrices
  atrisk <- outer(times, dtimes, ">=")
  hazinc <- as.vector(1 / (expci %*% atrisk))
  matrixP <- outer(expci, hazinc, "*") * atrisk
  matrixPO <- matrixP %*% t(matrixO)
  matrixM <- matrix(d, n, nd, byrow=FALSE) * (!atrisk) - matrixPO %*% (!atrisk)
  matrixW <- -crossprod(t(matrixPO))
  diag(matrixW) <- diag(matrixW) + rowSums(matrixPO)
  if (m > 0) {
    ZW <- crossprod(Z, matrixW)
    ZWZ <- ZW %*% Z
    ZWZinvZ <- solve(ZWZ, t(Z))
    X <- X - Z %*% solve(ZWZ, ZW) %*% X
    csm <- colSums(X*X)
    X[,csm < max(csm)*1e-14] <- 0
  }

  # martingale residuals
  Y <- rowSums(matrixPO) - d

  if (perms > 0) {
    # unlike the glm and lm models, we make permutations of the index vector rather than Y itself
    # This is for also being able to permute tr(XXW)
    piy <- match(unique(Y),Y)[match(Y, unique(Y))]  # sets indices of equal Y to be equal
    # all permutations if possible
    if (npermutations(piy) <= perms) {
      random <- FALSE
      permIY <- allpermutations(piy)
    } else {
      # otherwise random permutations
      random <- TRUE
      permIY <- replicate(perms, piy[sample(n)])
    }
    permY <- matrix(Y[permIY], nrow=length(Y))
  }

  # The test statistic
  if (perms == 0) {
    test <- function(subset, weights, calculateP = TRUE) {   # calculateP not used in survival models
      if (!missing(subset))
        X <- X[,subset,drop=FALSE]
      p <- ncol(X)
      if (p == 0) {
        return(c(p = NA, S = NA, ES = NA, sdS = NA, ncov=0))
      } else {
        if(!missing(weights))
          X <- X * matrix(ifelse(weights>0, sqrt(weights), -sqrt(weights)), n, p, byrow = TRUE)
  
        XX <- crossprod(t(X))
        if (dir) XX <- XX + dir * outer(rowSums(X), rowSums(X))
        Q <- crossprod(Y, XX) %*% Y
  
        # adjust Q in case of tied survival times
        if (ties) {
          tiecorrect <- sum( sapply(1:length(dtimes), function(i) {
            if (sum(matrixO[,i]) == 1)
              0
            else {
              matrixtie <- outer(matrixO[,i], matrixO[,i])
              diag(matrixtie) <- 0
              sum(XX[as.logical(matrixtie)]) - 2 * sum(matrixP[,i] %*% XX %*% matrixtie) + sum(matrixtie) * (matrixP[,i] %*% XX %*% matrixP[,i])
            }
          }))
          Q <- Q - tiecorrect
        }
  
        # calculate mean and variance of Q
        EQ <- sum(XX * matrixW)
        tussen <- matrix(diag(XX), n, nd) + 2 * XX %*% (matrixM - matrixP)
        matrixT <- tussen - matrix(colSums(matrixP * tussen), n, nd, byrow = TRUE)
        varQ <- sum(matrixPO * ((matrixT * matrixT) %*% t(matrixO)))
  
        # calculate p-value
        Z <- (Q - EQ) / sqrt(varQ)
        p.value <- pnorm(Z, lower.tail = FALSE)
        norm.const <- (n-m) * EQ / 100
  
        # positive = negative association with occurence of the event
  
        # give back
        return(c(p = p.value, S = Q/norm.const, ES = EQ/norm.const, sdS=sqrt(varQ)/norm.const, ncov=p))
      }
    }
  } else {
    test <- function(subset, weights, calculateP = TRUE) {     # calculateP not used in permutation testing
      if (!missing(subset))
        X <- X[,subset,drop=FALSE]
      p <- ncol(X)
      if (p == 0) {
        return(c(p = NA, S = NA, ES = NA, sdS = NA, ncov=0))
      } else {
        if(!missing(weights))
          X <- X * matrix(ifelse(weights>0, sqrt(weights), -sqrt(weights)), n, p, byrow = TRUE)
  
        XX <- crossprod(t(X))
        if (dir) XX <- XX + dir * outer(rowSums(X), rowSums(X))
        Q <- drop(crossprod(Y, XX) %*% Y)
        EQ <- sum(XX * matrixW)
  
        permQ <- colSums(permY * (XX %*% permY))
        permEQ <- apply(permIY, 2, function(prm) {
          pW <- matrixW[prm,prm]
          sum(XX * pW)
        })
  
        norm.const <- (n-m) * EQ / 100
  
        # calculate the p-value
        p.value <- mean(Q - EQ <= permQ - permEQ)
  
        # mean and variance of S for z-score
        ES <- EQ + mean(permQ - permEQ)
        VarS <- var(permQ - permEQ)
  
        # give back
        return(c(p = p.value, S = Q/norm.const, ES = ES/norm.const, sdS = sqrt(VarS)/norm.const, ncov = p))
      }
    }
  }

  # Get the adjusted X matrix
  getX <- function(subset, weights) {
    if (!missing(subset))
      X <- X[,subset,drop=FALSE]
    p <- ncol(X)
    if(!missing(weights))
      X <- X * matrix(ifelse(weights>0, sqrt(weights), -sqrt(weights)), n, p, byrow = TRUE)
    X
  }


  # Subjectsplot function
  subjectplot <- function(subset, weights, mirror) {
    if (!missing(subset))
      X <- X[,subset,drop=FALSE]
    p <- ncol(X)
    if(!missing(weights))
      X <- X * matrix(ifelse(weights>0, sqrt(weights), -sqrt(weights)), n, p, byrow = TRUE)

    XX <- crossprod(t(X))
    if (dir) XX <- XX + dir * outer(rowSums(X), rowSums(X))
    if (mirror) {
      R <- drop(XX %*% Y * sign(Y))
      ER <- diag(XX) * Y * sign(Y)
    } else {
      R <- XX %*% Y
      ER <- diag(XX) * Y
    }

    norm.const <- mean(abs(ER))

    XXmd <- XX
    diag(XXmd) <- 0
    if (m > 0) XXmd <- XXmd - XXmd %*% matrixW %*% Z %*% ZWZinvZ
    varR <- diag(XXmd %*% matrixW %*% t(XXmd))

    cbind(R=R/norm.const, ER=ER/norm.const, sdR=sqrt(varR)/norm.const, Y=Y)
  }

  # permutations for permutations histogram
  permutations <- function(subset, weights) {
      if (!perms)
        stop("histogram only for permutation version of the test")

      if (!missing(subset))
        X <- X[,subset,drop=FALSE]
      p <- ncol(X)
      if(!missing(weights))
        X <- X * matrix(ifelse(weights>0, sqrt(weights), -sqrt(weights)), n, p, byrow = TRUE)

      XX <- crossprod(t(X))
      if (dir) XX <- XX + dir * outer(rowSums(X), rowSums(X))
      Q <- drop(crossprod(Y, XX) %*% Y)
      EQ <- sum(XX * matrixW)

      permQ <- colSums(permY * (XX %*% permY))
      permEQ <- apply(permIY, 2, function(prm) {
        pW <- matrixW[prm,prm]
        sum(XX * pW)
      })

      norm.const <- (n-m) * EQ / 100

      # calculate the p-value
      permS <- EQ + permQ - permEQ

      return(list(S = Q/norm.const, permS = permS/norm.const))
  }

  df <- function() {
    c(n,m,ncol(X))
  }

  nperms <- function() {
    if (perms)
      c(n = ncol(permY), random = random)
    else
      c(0, FALSE)
  }

  cov.names <- function(subset)
    if (missing(subset))
      colnames(X)
    else
      colnames(X[,subset,drop=FALSE])

  dist.cor <- function(subset, weights, transpose = FALSE) {
    if (!missing(subset))
      X <- X[,subset,drop=FALSE]
    p <- ncol(X)
    if(!missing(weights))
      X <- X * matrix(ifelse(weights>0, sqrt(weights), -sqrt(weights)), n, p, byrow = TRUE)
    if (transpose) 
      X <- t(X)
    sds <- sqrt(colSums(X*X))
    sds[sds < max(sds)*1e-14] <- Inf
    crossprod(X) / outer(sds,sds)
  }
  
  positive <- function(subset, weights) {
    if (!missing(subset))
      X <- X[,subset,drop=FALSE] 
    xy <- crossprod(X, Y)
    drop(xy >= 0)+0
  }
  
  getY <- function() Y

  return(list(test=test,
    getX = getX,
    cov.names = cov.names,
    cor=dist.cor,
    getY = getY,
    subjectplot = subjectplot,
    permutations = permutations,
    nperms = nperms,
    positive = positive,
    df=df))

}


##############################################
# Multinomial logistic model
##############################################
.multinomialtest <- function(response, Z, X, dir, perms) {

  # establish dimensions of the problem
  m <- ncol(Z)
  n <- nrow(Z)

  # fit the null model
  form <- formula(.makeFormula(m>0, FALSE))
  null.fit <- mlogit(form)
  fitted <- fitted.values(null.fit)
  G <- ncol(fitted)
  if ((perms > 0) && (!all(sapply(1:G, function(i) all(fitted[,G] == fitted[1,G])))))
    stop("Multinomial model with covariates:\n\tPermutation testing is not possible.", call.=FALSE)
  bY <- as.vector(null.fit@y - fitted)
  cY <- apply(null.fit@y==1, 1, which)
  bX <- diag(G) %x% X
  W <- matrix(0,n*G,n*G)
  for (i in 1:G) 
    for (j in 1:G)
      diag(W[(i-1)*n + 1:n, (j-1)*n + 1:n]) <- if (i==j) fitted[,i]*(1-fitted[,i]) else -fitted[,i]*fitted[,j]

  if (perms == 0) {
    eW <- eigen(W,symmetric=TRUE)
    eW$vectors <- eW$vectors[,1:(n*(G-1)),drop=FALSE]
    eW$values <- eW$values[1:(n*(G-1))]
    halfW <- eW$vectors %*% diag(eW$values^(1/2)) %*% t(eW$vectors)   
    # adjust the alternative
    if (m > 0) {
      bZ <- diag(G) %x% Z       
      bZW <- crossprod(bZ, W)
      eZWZ <- eigen(bZW %*% bZ, symmetric=TRUE)    
      eZWZ$vectors <- eZWZ$vectors[,1:(m*(G-1)),drop=FALSE]
      eZWZ$values <- eZWZ$values[1:(m*(G-1))]
      ZWZinvZW <- eZWZ$vectors %*% diag(1/eZWZ$values, nrow=m*(G-1)) %*% crossprod(eZWZ$vectors, bZW) 
      bX <- bX - bZ %*% (ZWZinvZW %*% bX) 
      csm <- colSums(bX*bX)
      bX[,csm < max(csm)*1e-14] <- 0
    }
  } else {
    equal.size.groups <- all(table(cY)==table(cY)[1])
    if (equal.size.groups)
      nperms <- npermutations(cY[-1])
    else
      nperms <- npermutations(cY)
    # all permutations if possible
    if (nperms <= perms) {
      random <- FALSE
      if (equal.size.groups)
        permY <- rbind(cY[1], allpermutations(cY[-1]))
      else
        permY <- allpermutations(cY)
    } else {         
      # otherwise random permutations
      random <- TRUE
      permY <- replicate(perms, cY[sample(n)])
    }
    permY <- apply(permY, 2, function(prm) {
      as.vector(sapply(1:G, function(i) { prm == i }) - fitted)
    }) 
    # adjust the alternative
    if (m > 0) {
      bZ <- diag(G) %x% Z       
      bZW <- crossprod(bZ, W)
      eZWZ <- eigen(bZW %*% bZ, symmetric=TRUE)    
      eZWZ$vectors <- eZWZ$vectors[,1:(m*(G-1)),drop=FALSE]
      eZWZ$values <- eZWZ$values[1:(m*(G-1))]
      ZWZinvZW <- eZWZ$vectors %*% diag(1/eZWZ$values,nrow=m*(G-1)) %*% crossprod(eZWZ$vectors, bZW) 
      bX <- bX - bZ %*% (ZWZinvZW %*% bX) 
      csm <- colSums(bX*bX)
      bX[,csm < max(csm)*1e-14] <- 0
    }
  }

  # utility to select columns of bX
  select.help <- 1:ncol(X)
  names(select.help) <- colnames(X)
  selector <- function(set) outer(ncol(X)*(1:G-1), select.help[set], "+")

  norm.const <- (n-m) / 100

  # calculate the test statistic
  # 1: asymptotic version
  if (perms == 0) {
    test <- function(subset, weights, calculateP = TRUE) {
      if (!missing(subset))
        bX <- bX[,selector(subset),drop=FALSE]
      p <- ncol(bX)/G
      if (p == 0) {
        return(c(p = NA, S = NA, ES = NA, sdS = NA, ncov=0))
      } else {
        if(!missing(weights))
          bX <- bX * matrix(rep(ifelse(weights>0, sqrt(weights), -sqrt(weights)),G), n*G, p*G, byrow = TRUE)
  
        XX <- crossprod(t(bX))
        if (dir) 
          for (i in 1:G) { 
            rs <- rowSums(bX[,p*(i-1)+1:p])
            XX <- XX + dir * outer(rs, rs)
          }
        dXX <- XX
        dXX[(row(XX) - col(XX)) %% n != 0] <- 0
                     
        S <- sum(bY * (XX %*% bY)) / sum(bY * (dXX %*% bY))
  
        # get lambdas
        if (m > 0) {
          dXX <- dXX - bZ %*% (ZWZinvZW %*% dXX)
          dXX <- dXX - dXX %*% t(ZWZinvZW) %*% t(bZ)
        }
        XX <- halfW %*% XX %*% halfW
        dXX <- halfW %*% dXX %*% halfW
        if (is.nan(S))
          lams <- 0
        else
          lams <- eigen(XX - S*dXX, symmetric=TRUE, only.values=TRUE)$values
  
        # calculate the p-value
        if (calculateP)
          p.value <- .getP(lams)
        else
          p.value <- NA
  
        # mean and variance of S for z-score
        # use series approximation as given in Paolella (2003) 319:
        mu.num <- sum(diag(XX))
        mu.den <- sum(diag(dXX))
        var.num <- 2*sum(XX*XX)
        var.den <- 2*sum(dXX*dXX)
        cov.term <- 2*sum(XX*dXX)
  
        ES <- (mu.num/mu.den) * (1 - cov.term/(mu.num*mu.den) + var.den/(mu.den^2))
        VarS <- (mu.num^2/mu.den^2) * (var.num/(mu.num^2)  + var.den/(mu.den^2) - 2*cov.term/(mu.num*mu.den))
  
        # give back
        return(c(p = p.value, S = S/norm.const, ES = ES/norm.const, sdS = sqrt(VarS)/norm.const, ncov=p))
      }
    }
  } else {
  # 2: permutation version
    test <- function(subset, weights, calculateP=TRUE) {  # calculateP not used in permutation testing
      if (!missing(subset))
        bX <- bX[,selector(subset),drop=FALSE]
      p <- ncol(bX)/G
      if (p == 0) {
        return(c(p = NA, S = NA, ES = NA, sdS = NA, ncov=0))
      } else {
        if(!missing(weights))
          bX <- bX * matrix(rep(ifelse(weights>0, sqrt(weights), -sqrt(weights)),G), n*G, p*G, byrow = TRUE)
          
        XX <- crossprod(t(bX))
        if (dir) 
          for (i in 1:G) { 
            rs <- rowSums(bX[,p*(i-1)+1:p])
            XX <- XX + dir * outer(rs, rs)
          }
        dXX <- XX
        dXX[(row(XX) - col(XX)) %% n != 0] <- 0
      
        S <- sum(bY * (XX %*% bY)) / sum(bY * (dXX %*% bY))
        permS <- colSums(permY * (XX %*% permY)) / colSums(permY * (dXX %*% permY))
  
        # mean and variance of S for z-score
        ES <- mean(permS)
        VarS <- var(permS)
  
        # calculate the p-value
        p.value <- mean(S <= permS)
  
        # give back
        return(c(p = p.value, S = S/norm.const, ES = ES/norm.const, sdS = sqrt(VarS)/norm.const, ncov=p))
      }
    }
  }

  # Get the adjusted X matrix
  getX <- function(subset, weights) {
    if (!missing(subset))
      bX <- bX[,selector(subset),drop=FALSE]
    p <- ncol(bX)/G
    if(!missing(weights))
      bX <- bX * matrix(rep(ifelse(weights>0, sqrt(weights), -sqrt(weights)),G), n*G, p*G, byrow = TRUE)
    bX
  }

  # Subjectsplot function
  subjectplot <- function(subset, weights, mirror) {
    if (!missing(subset))
      bX <- bX[,selector(subset),drop=FALSE]
    p <- ncol(bX)/G
    if(!missing(weights))
      bX <- bX * matrix(rep(ifelse(weights>0, sqrt(weights), -sqrt(weights)),G), n*G, p*G, byrow = TRUE)

    XX <- crossprod(t(bX))
    if (dir) 
      for (i in 1:G) { 
        rs <- rowSums(bX[,p*(i-1)+1:p])
        XX <- XX + dir * outer(rs, rs)
      }    
    dXX <- XX
    dXX[(row(XX) - col(XX)) %% n != 0] <- 0

    mY <- matrix(bY, n*G, n)
    mY[(row(mY) - col(mY)) %% n != 0] <- 0
    if (mirror) {
      R <- drop(crossprod(mY, XX %*% bY))
      ER <- drop(crossprod(mY, dXX %*% bY))
    } else 
      stop("mirror=FALSE not appropriate for the multinomial model")

    norm.const <- mean(abs(ER))

    XXmd <- XX - dXX
    if (m > 0) 
      XXmd <- XXmd - XXmd %*% t(ZWZinvZW) %*% t(bZ)
    XXmdW <- crossprod(mY, XXmd) %*% halfW 
    varR <- drop(rowSums(XXmdW * XXmdW))

    cbind(R=R/norm.const, ER=ER/norm.const, sdR=sqrt(varR)/norm.const, Y=cY)
  }

  # permutations for permutations histogram
  permutations <- function(subset, weights) {
      if (!perms)
        stop("histogram only for permutation version of the test")

      if (!missing(subset))
        bX <- bX[,selector(subset),drop=FALSE]
      p <- ncol(bX)/G
      if(!missing(weights))
        bX <- bX * matrix(rep(ifelse(weights>0, sqrt(weights), -sqrt(weights)),G), n*G, p*G, byrow = TRUE)
        
      XX <- crossprod(t(bX))
      if (dir) 
        for (i in 1:G) { 
          rs <- rowSums(bX[,p*(i-1)+1:p])
          XX <- XX + dir * outer(rs, rs)
        }
      dXX <- XX
      dXX[(row(XX) - col(XX)) %% n != 0] <- 0
    
      S <- sum(bY * (XX %*% bY)) / sum(bY * (dXX %*% bY))
      permS <- colSums(permY * (XX %*% permY)) / colSums(permY * (dXX %*% permY))

      return(list(S = S/norm.const, permS = permS/norm.const))
  }

  df <- function() {
    c(n,m,ncol(X))
  }

  nperms <- function() {
    if (perms)
      c(n = ncol(permY), random = random)
    else
      c(0, FALSE)
  }

  cov.names <- function(subset)
    if (missing(subset))
      colnames(X)
    else
      colnames(X[,subset,drop=FALSE])

  dist.cor <- function(subset, weights, transpose=FALSE) {
    if (!missing(subset))
      X <- X[,subset,drop=FALSE]
    p <- ncol(X)
    if(!missing(weights))
      X <- X * matrix(ifelse(weights>0, sqrt(weights), -sqrt(weights)), n, p, byrow = TRUE)
    if (transpose)
      uit <- matrix(0,n,n)
    else
      uit <- matrix(0,p,p)
    for (i in 1:G) {
      rngN <- 1:n + (i-1)*n
      if (m > 0) {
        ZW <- Z * matrix(diag(W[rngN,rngN]), nrow(Z), ncol(Z))
        adjX <- X - Z %*% solve(crossprod(ZW, Z), crossprod(ZW, X))
      } else
        adjX <- X
      if (transpose)
        adjX <- t(adjX)
      sds <- sqrt(colSums(adjX*adjX))
      sds[sds < max(sds)*1e-14] <- Inf
      uit <- uit + crossprod(adjX) / outer(sds,sds)
    }
    uit/G
  }

  
  positive <- function(subset) {
    if (!missing(subset))
      X <- X[,subset,drop=FALSE] 
    drop(apply(cor(X, null.fit@y - fitted), 1, which.max))
  }  
          

  
  getY <- function() bY

  return(list(test=test,
    getX = getX,
    cov.names = cov.names,
    cor=dist.cor,
    getY = getY,
    subjectplot = subjectplot,
    permutations = permutations,
    nperms = nperms,
    positive = positive,
    df=df))


}
