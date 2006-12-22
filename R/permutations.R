#==========================================================
# Histogram method for permutations of "gt.result"
#==========================================================
setMethod("hist", "gt.result",
            function(x, ...) 
{
  gt.hist <- function(x, geneset = NULL, ...) {
    if (!is.null(geneset))
      x <- x[geneset]
    if (.nPathways(x) > 1)
      stop("more than one geneset in x", call. = FALSE)
    if (ncol(x@PermQs) == 0)
      stop("no permutations in x: try hist(permutations(x)) instead", call. = FALSE)
      
    Qs <- as.vector(x@PermQs)
    Q <- .Q(x)
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
  if (!(.model(gt) == "survival") && .adjusted(gt))
    stop("the permutation procedure is not yet implemented for the adjusted global test", call. = FALSE)      
      
  # check correct input of nperm
  if ( !(nperm > 0) )
    stop("option nperm should be a positive integer", call. = FALSE)

  if (min(nperm, .nPerms(gt)) > ncol(gt@PermQs)) {
    # extract the right test.genes vector
    if (!missing(geneset))
      gt <- gt[geneset]
    nTested <- .nTested(gt)
    realQ <- .Q(gt)
    n <- .nSamples(gt)
    model <- .model(gt)
    adjusted <- .adjusted(gt)
  
    # recreate Y and necessary matrices
    if (model == "survival") {
      times <- as.vector(fit(gt)$y)[1:n]
      d <- as.vector(fit(gt)$y)[(n+1):(2*n)] 
      expci <- exp(fit(gt)$linear.predictors)
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
    } else if (model == "linear") {
      Y <- fit(gt)$y - fitted.values(fit(gt))
      mu2 <- sum(Y*Y)/df.residual(fit(gt))
    } else if (model == "logistic") {
      Y <- fit(gt)$y - fitted.values(fit(gt))
      mu2 <- fit(gt)$weights[1]
    } else if (model == "multinomial") {
      Y <- .Y(gt)
    }
  
    # Make the permutations of Y
    if (nperm >= .nPerms(gt)) { # Use all possible permutations
      gt@method <- 4
      gt@PermQs <- matrix(,.nPathways(gt),0)
      if ((model == "logistic") && !adjusted) {
        m <- sum(Y == Y[1])
        if (m == n/2) {
          Y.pm <- rbind(TRUE, .allperms2(m-1,n-1))
        } else {
          Y.pm <- .allperms2(m,n)
        }
        Y.pm <- Y.pm - m/n
      } else if (model == "linear" && !adjusted) {
        Y.pm <- apply(.allperms(1:n), 2, function(pm) Y[pm])
      } else if (model == "multinomial") {
        Y <- apply(.Y(gt), 1, which.max)
        counts <- sapply(unique(Y), function(ix) sum(Y == ix))
        g <- length(counts)
        pms <- .allpermsG(counts, counts)
        Y.pm <- lapply(as.list(1:g), function(ix) {
          temp <- (pms == ix)
          temp - colMeans(temp)
        })
      } else if (model == "survival") {
        if (!ties) {
          Y.pm <- apply(.allperms(1:n), 2, function(pm) Y[pm])
        } else {
          pms <- .allperms(1:n)
        }
      }
    } else { # Use random permutations
      gt@method <- 5
      nperm <- nperm - ncol(gt@PermQs)
      pms <- apply( matrix(rnorm(n * nperm), n, nperm), 2, sort.list )
      if (model == "survival") {
        if (!ties) {
          Y.pm <-  apply(pms, 2, function(pm) Y[pm])
        }
      } else if (model == "multinomial") {
        Y <- lapply(as.list(1:ncol(Y)), function(ix) Y[,ix])
        Y.pm <- lapply(Y, function(yy) {
          apply(pms, 2, function(pm) yy[pm])
        })
      } else {
        Y.pm <- apply(pms, 2, function(pm) Y[pm])
      }
    }
    
    QQs <- t(sapply(1:.nPathways(gt), function(index) {
      if (nTested[index] == 0)
        Qs <- rep(0, times = nperm)
      else {
        # Recreate XX
        X <- gt@eX[gt@genesets[[index]], , drop = FALSE]
        m <- nrow(X)
        XX <- crossprod(X) / m
        Q <- realQ[index]
    
        # Calculate Q for nperm permutations of Q
        if (model %in% c("logistic", "linear")) {
          Qs <- colSums(( XX %*% Y.pm ) * Y.pm) / mu2
        } else if (model == "multinomial") {
          Qs <- rowSums(sapply(Y.pm, function(yy) colSums(( XX %*% yy ) * yy)))
        } else { 
          # survival model:
          if (adjusted) {
            R <- crossprod(.IminH(gt), XX) %*% .IminH(gt)
          } else {
            R <- XX
          }
          if (!ties) {
            Q.pms <- colSums(( XX %*% Y.pm ) * Y.pm)
            EQs <- rep(0, nperm)
            varQs <- rep(0, nperm)
            for (j in 1:nd) {
              pj <- matrixP[,j]
              PJ <- matrix( pj[pms], n, nperm )
              mj <- matrixM[,j]
              MJ <- matrix( mj[pms], n, nperm )
              tussen <- matrix(diag(R), n, nperm) + 2 * R %*% (MJ - PJ)
              TJ <- tussen - matrix(colSums(PJ * tussen), n, nperm, byrow = TRUE)
              varQs <- varQs + colSums(PJ * TJ * TJ)
              EQs <- EQs + diag(R) %*% PJ - colSums(( R %*% PJ ) * PJ)
            }
            Qs <- (Q.pms - EQs) / sqrt(varQs)
          } else {
            # ties present: much slower generation of permutations
            Qs <- sapply(1:nperm, function(ix) {
              R <- R[pms[,ix], pms[,ix]]
              XX <- XX[pms[,ix], pms[,ix]]
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
    gt@PermQs <- cbind(gt@PermQs, QQs)
    EQ <- apply(gt@PermQs, 1, mean)
    seQ <- apply(gt@PermQs, 1, sd)
    p.value <- apply(gt@PermQs > .Q(gt) * 0.999999, 1, mean)
    gt@res[,4:6] <- cbind(EQ, seQ, p.value)
    gt@res <- gt@res[,1:6, drop = FALSE]
  }
  gt
}
#==========================================================


#==========================================================
# Lists all permutations for the two-group case
#==========================================================
.allperms2 <- function(m, n) {
  if (n == 1) {
    if (m == 0) {
      app <- FALSE
    } else if (m==1) {
      app <- TRUE
    }
  } else {
    total <- choose(n,m)
    top <- choose(n-1,m-1)
    bottom <- choose(n-1,m)
    app <- matrix(,n,choose(n,m))
    app[1,] <- c(rep(TRUE, top), rep(FALSE, bottom))
    if (m > 0)
      app[2:n,1:top] <- .allperms2(m-1,n-1)
    if (m < n)
      app[2:n,(top+1):total] <- .allperms2(m,n-1)
  }  
  app
}


#==========================================================
# Lists all permutations for the multiple-group case
#==========================================================
.allpermsG <- function(counts, grouping) {
  n <- sum(counts)
  if (n == 1) {
    app <- which.max(counts)
  } else {
    total <- .nPermsG(counts, grouping)
    app <- matrix(,n,total)
    choosable <- (counts > 0) & (is.na(grouping) | (1:length(counts) %in% match(unique(grouping[!is.na(grouping)]), grouping)))
    choosable <- (1:length(counts))[choosable]
    ix <- 0
    for (iy in choosable) {
      countstemp <- counts
      countstemp[iy] <- counts[iy] - 1
      groupingtemp <- grouping
      groupingtemp[iy] <- NA
      size <- .nPermsG(countstemp, groupingtemp)
      app[1,(ix+1):(ix+size)] <- iy
      app[2:n, (ix+1):(ix+size)] <- .allpermsG(countstemp, groupingtemp)
      ix <- ix + size
    }
  }  
  app
}

.nPermsG <- function(counts, grouping) {
  total <- .mchoose(counts)
  if (any(!is.na(grouping))) {
    correction <- prod(factorial(sapply(unique(grouping[!is.na(grouping)]), function(cc) sum(grouping == cc, na.rm=TRUE))))
  } else {
    correction <- 1
  }
  total / correction
}

.mchoose <- function(counts) {
  out <- choose(sum(counts), counts[1])
  if (length(counts) > 2)
    out <- out * .mchoose(counts[-1])
  out
}


#==========================================================
# Lists all permutations for the general case
#==========================================================
.allperms <- function(nums) {
  n <- length(nums)
  if (n == 1) {
    app <- nums
  } else {
    app <- matrix(,n,factorial(n))
    for (ix in 1:length(nums)) {
      range <- 1:factorial(n-1) + (ix - 1) * factorial(n-1)
      app[1,range] <- nums[ix]
      app[2:n,range] <- .allperms(nums[-ix])
    }
  }
  app
}


#==========================================================
# Calculates the number of permutations
#==========================================================
.nPerms <- function(gt) {
  n <- .nSamples(gt)
  if ((.model(gt) == "logistic") && !.adjusted(gt)) {
    Y <- .Y(gt)
    m <- sum(Y == Y[1])
    if (m == n/2) {
      out <- choose(n,m) / 2
    } else {
      out <- choose(n,m)
    }
  } else if ((.model(gt) == "multinomial") && !.adjusted(gt)) {
    Y <- apply(.Y(gt), 1, which.max)
    counts <- sapply(unique(Y), function(ix) sum(Y == ix))
    out <- .nPermsG(counts, counts)
  } else {
    out <- ifelse(n<=100, factorial(n), Inf)
  }
  out
}
