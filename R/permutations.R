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
    stop("the permutation procedure is not applicable for the adjusted global test", call. = FALSE)      
      
  # check correct input of nperm
  if ( !(nperm > 0) )
    stop("option nperm should be a positive integer", call. = FALSE)
  if (nperm <= ncol(gt@PermQs))
    gt@PermQs <-  gt@PermQs[,1:nperm]
  else {
    nperm <- nperm - ncol(gt@PermQs)
    # extract the right test.genes vector
    if (!missing(geneset))
      gt <- gt[geneset]
    nTested <- .nTested(gt)
    realQ <- .Q(gt)
    n <- .nSamples(gt)
    model <- .model(gt)
    adjusted <- .adjusted(gt)
    
    QQs <- t(sapply(1:.nPathways(gt), function(index) {
      if (nTested[index] == 0)
        Qs <- rep(0, times = nperm)
      else {
        # Recreate XX
        X <- gt@eX[gt@genesets[[index]], , drop = FALSE]
        m <- nrow(X)
        XX <- crossprod(X) / m
        Q <- realQ[index]

        # recreate Y and necessary matrices
        if (model == "survival") {
          if (adjusted) {
            R <- crossprod(.IminH(gt), XX) %*% .IminH(gt)
          } else {
            R <- XX
          }
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
        }        

        # Calculate Q for nperm permutations of Q
        if (model != 'survival') {
          # Recalculate the Q-value for permutations of Y
          permQ <- function(npm) {
            pms <- apply( matrix(rnorm(n * npm), n, npm), 2, sort.list )
            Y.pm <-  matrix( Y[pms], n, npm )
            colSums(( XX %*% Y.pm ) * Y.pm) / mu2
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
    gt@PermQs <- cbind(gt@PermQs, QQs)
  }
  EQ <- apply(gt@PermQs, 1, mean)
  seQ <- apply(gt@PermQs, 1, sd)
  p.value <- apply(gt@PermQs >= matrix(.Q(gt), nrow(gt@PermQs), ncol(gt@PermQs)), 1, mean)
  gt@res[,4:6] <- cbind(EQ, seQ, p.value)
  gt
}
#==========================================================
