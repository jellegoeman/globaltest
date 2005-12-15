#==========================================================
# Sampleplot plots the influence of each sample on the outcome 
#   of the test statistic
# See help(sampleplot) for details
#==========================================================
sampleplot <- function(gt, geneset, samplesubset, plot = TRUE, scale = TRUE, drawlabels = TRUE, labelsize = 0.6,...)
{
  # check correct input of gt
  if ( !is(gt, "gt.result"))
      stop("geneplot should be applied to a globaltest result", call. = FALSE)
      
  # extract the right test.genes vector
  if (!missing(geneset))
    gt <- gt[geneset]
  if (length(gt) > 1)
    stop("more than one geneset selected", call. =FALSE)
  if (.nTested(gt) == 0)
    stop("empty pathway", call. = FALSE)
  geneset <- gt@genesets[[1]]

  # reconstruct X and Y
  X <- as.matrix(gt@eX[geneset,])
  n <- ncol(X)
  if (missing(samplesubset))
    samplesubset <- 1:n
  else {  
    if (all(samplesubset %in% rownames(X)))
      samplesubset <- match(samplesubset, rownames(X))
    else
      if (!all(samplesubset %in% 1:n) & !all(samplesubset %in% rownames(X)))
        stop("samplesubset should contain names or numbers of samples", call. = FALSE)
  }
  m <- nrow(X)
  model <- .model(gt)
  adjusted <- .adjusted(gt)
  if (adjusted) {
    adjX <- X %*% .IminH(gt)
  } else {
    adjX <- X
  }
  
  # calculate influence per sample and expected influence
  if (model != 'survival') {
    Y <- .Y(gt)
    if ((model == "logistic") && !adjusted) {
      mu2 <- fit(gt)$weights[1]
      Y <- Y / sqrt(mu2)
    } else if (model == "linear") {
      mu2 <- sum(Y*Y)/df.residual(fit(gt))
      Y <- Y / sqrt(mu2)
    }  
    if (adjusted) {
      XXY <- crossprod(.IminH(gt), crossprod(X) %*% Y)
    } else {
      XXY <- crossprod(X) %*% Y
    }
    influence <- (n/m) * Y * XXY
    up <- (sign(Y) == 1) 
    R <- crossprod(adjX)
    if (model == 'logistic' && adjusted) {
      mu2 <- fit(gt)$weights
      RR <- (R * R) %*% diag(mu2)
    } else {
      RR <- R * R
    }
    Einf <- (n/m) * Y * Y * colSums(adjX*adjX)
    RR <- matrix(t(RR)[diag(n) == 0], n, n-1, byrow = TRUE)
    varinf <- (n/m)^2 * Y * Y * rowSums(RR)
    sd.inf <- sqrt(varinf)
  } else {
    # survival model
    R <- crossprod(adjX)
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
    matrixM <- diag(d) %*% outer(times, dtimes, "<") - matrixPO %*% outer(times, dtimes, "<")
    matrixW <- diag(rowSums(matrixPO)) - matrixPO %*% t(matrixPO)
    Y <- rowSums(matrixPO) - d
    influence <- colSums(crossprod(adjX, X) * (outer(Y,Y) - matrixW))
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
          XX <- crossprod(X) / n
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
  colnames(res) <- c("Influence", "Expected", "SD", "UP")
  rownames(res) <- colnames(X)
  if (model == "linear") {
    nameY <- as.character(.formula(gt)[[2]])
    colourCode <- c("+", "-")
    if (adjusted) {
      legend <- c(paste("positive residual", nameY),
        paste("negative residual", nameY))
    } else {
      legend <- c(paste("large values of", nameY),
        paste("small values of", nameY))
    }
  } else if (model == 'logistic') {
    colourCode <- c(.levels(gt)[2], .levels(gt)[1] )
    legend <- c(paste(.levels(gt)[2], "samples"), paste(.levels(gt)[1], "samples"))
  } else if (model == 'survival') {
    colourCode <- c("late", "early")
    legend <- c("late event time or censored", "early event time")
  }


  gtbar <- new("gt.barplot",     
    res = res,
    labelsize = labelsize, 
    drawlabels = drawlabels,
    colourCode = colourCode,
    legend = legend)
  gtbar <- gtbar[samplesubset]
  if (scale)
    gtbar <- scale(gtbar)
  if (plot) {
    plot(gtbar, ...)
    invisible(gtbar)
  } else {
    gtbar
  }
}
#==========================================================
