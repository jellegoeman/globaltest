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
  X <- as.matrix(gt@eX[geneset,,drop=FALSE])
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
  
  # calculate influence per sample and expected influence
  if (model %in% c("linear", "logistic")) {
    Y <- .Y(gt)
    if ((model == "logistic") && !adjusted) {
      mu2 <- fit(gt)$weights[1]
      Y <- Y / sqrt(mu2)
    } else if (model == "linear") {
      mu2 <- sum(Y*Y)/df.residual(fit(gt))
      Y <- Y / sqrt(mu2)
    }  
    if (model == "logistic") {
       Y <- -Y
    }
    if (adjusted) {
      adjX <- X %*% .IminH(gt)
      XXY <- crossprod(.IminH(gt), crossprod(X) %*% Y)
    } else {
      adjX <- X
      XXY <- crossprod(X) %*% Y
    }
    influence <- as.vector((n/m) * Y * XXY)
    up <- 1 + (sign(Y) == 1) 
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
  } else if ((model == "multinomial") && !adjusted) {
    Y <- .Y(gt)
    g <- ncol(Y)
    mu <- fitted.values(fit(gt))[1,]
    matrixW <- -mu %o% mu + diag(mu)
    XX <- crossprod(X) / m
    influence <- n * rowSums(sapply(1:g, function(out) Y[,out] * (XX %*% Y[,out])))
    Einf <- n * rowSums(Y * Y * (diag(XX) %o% rep(1,g)))
    RR <- matrix(t(XX * XX)[diag(n) == 0], n, n-1, byrow = TRUE)
    varinf <- n^2 * rowSums(Y * (Y %*% matrixW)) * rowSums(RR)
    sd.inf <- sqrt(varinf)
    up <- apply(Y, 1, which.max)
  } else if ((model == "multinomial") && adjusted) {
    Y <- .Y(gt)
    bigY <- as.vector(Y)
    g <- ncol(Y)
    range <- lapply(as.list(1:g), function(ix) ((ix-1)*n+1):(ix*n))
    mu <- fitted.values(fit(gt))
    matrixW <- matrix(0,g,n*g)
    for (ix in 1:g) 
      matrixW[,range[[ix]]] <- -t(mu)
    matrixW <- matrixW * (rep(1,g) %o% as.vector(mu))
    for (ix in 1:g) 
      matrixW[ix,range[[ix]]] <- matrixW[ix,range[[ix]]] + mu[,ix]
    bigXX <- diag(g) %x% crossprod(X) 
    R <- crossprod(.IminH(gt), bigXX) %*% .IminH(gt)
    XXY <- matrix(crossprod(.IminH(gt), bigXX %*% bigY), n, g)
    influence <- (n/m) * rowSums(Y * XXY)
    Einf <- (n/m) * rowSums(sapply(1:g, function(out1) 
      rowSums(sapply(1:g, function(out2) Y[,out1] * Y[,out2] * diag(R[range[[out1]],range[[out2]]])))
    ))
    M <- sapply(1:n, function(j)
      rowSums(sapply(1:g, function(u)
        rowSums(sapply(1:g, function(v) 
          rowSums(Y * matrix(R[(u-1)*n+j,],n,g)) * rowSums(Y * matrix(R[(v-1)*n+j,],n,g)) * matrixW[u, (v-1)*n+j]
        ))
      ))
    )
    diag(M) <- 0
    varinf <- (n/m)^2 * rowSums(M)
    sd.inf <- sqrt(varinf)
    up <- apply(Y, 1, which.max)
  } else if (model == "survival") {
    # survival model
    if (adjusted) {
      adjX <- X %*% .IminH(gt)
    } else {
      adjX <- X
    }
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
    up <- 1 + (sign(Y) == 1) 
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
    colour <- c(2,3)
    colourCode <- c("-", "+")
    if (adjusted) {
      legend <- c(paste("negative residual", nameY), paste("positive residual", nameY))
    } else {
      legend <- c(paste("small values of", nameY), paste("large values of", nameY))
    }
  } else if (model == 'logistic') {
    colour <- c(2,3)
    colourCode <- .levels(gt)
    legend <- c(paste(.levels(gt)[1], "samples"), paste(.levels(gt)[2], "samples"))
  } else if (model == "multinomial") {
    colourCode <- .levels(gt)
    colour <- 1 + 1:length(.levels(gt))
    legend <- .levels(gt)
  } else if (model == 'survival') {
    colourCode <- c("early", "late")
    colour <- c(2,3)
    legend <- c("early event time", "late event time or censored")
  }


  gtbar <- new("gt.barplot",     
    res = res,
    labelsize = labelsize, 
    drawlabels = drawlabels,
    colour = colour,
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
