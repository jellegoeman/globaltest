.getP <- function(lams) {
  if (all(lams==0))
    p.value <- NaN
  else {
    # first try Imhof (quick but less accurate)
    try.p <- .pImhof(lams)
    if (is.na(try.p$value) || (try.p$value <= 10*try.p$error))  {
      p.value <- .genF(1,lams)
    } else
      p.value <- try.p$value
  }
  p.value
}


.weed <- function(lam, acc = 50)
{
  thresh <- 1/acc
  lam <- sort(lam,decreasing=TRUE)
  mmm <- length(lam)
  while ((mmm > 0) && (lam[mmm]/lam[1] < thresh)) {
    qqq <- mmm-1
    rrr <- mmm-2
    lam[qqq]  <-  lam[qqq] + lam[mmm]
    while ((rrr>0) && (lam[rrr]<lam[qqq])) {
      lam[rrr:qqq] <- mean(lam[rrr:qqq])
      rrr <- rrr-1
    }
    mmm <- qqq
  }
  lam <- lam[seq_len(mmm)]
  return(lam)
}

.getC <- function(lams, beta, eps = 1e-10) {
                                       
  lams <- sort(lams)
  ready <- FALSE
  ix <- 1
  d <- NULL
  c <- prod(sqrt(beta/lams))
  rest.c <- 1-c
  d.base <- 1-beta/lams
  m <- length(lams)
               
  while (!ready)
  {
    d <- c(d, 0.5*sum(d.base^ix))
    c <- c(c,mean(c[1:ix]*d[ix:1]))
    if (rest.c > 100*.Machine$double.neg.eps)
      rest.c <- rest.c-c[ix+1]
    else
      rest.c <- rest.c*lams[1]/lams[m]
    ready <- (rest.c<eps) || (c[ix+1] < 0)
    ix <- ix+1
  }
  if (rest.c > 0) 
    c <- c(c, rest.c)
  if (any(c<0))
    return(NULL)
  else
    return(c)
}

.ruben <- function(lams) {
  mx <- max(lams)
  mn <- min(lams)
  return(2 * mx * mn / (mx + mn))
}

.genF <- function(x, lams, eps = 1e-10, acc = c(50,50)) {

  lams.pos <- .weed(lams[lams>0], acc = acc[1])
  while (prod(sqrt(.ruben(lams.pos)/lams.pos)) < .Machine$double.xmin) {     # prevents c=Inf/Inf, sacrificing accuracy
    acc[1] <- mean(c(acc[1],1))
    lams.pos <- .weed(lams.pos, acc=acc[1])
  }
  m.pos <- length(lams.pos)
  lams.neg <- .weed(-lams[lams<0], acc = acc[2])  
  while (prod(sqrt(.ruben(lams.neg)/lams.neg)) < .Machine$double.xmin) {     # prevents c=Inf/Inf, sacrificing accuracy
    acc[2] <- mean(c(acc[2],1))
    lams.neg <- .weed(lams.neg, acc=acc[2])
    print(acc)
  }
  m.neg <- length(lams.neg)
                      
  if (m.neg == 0)
    p.value <- as.numeric(x>0)
  else if (m.pos == 0)
    p.value <- as.numeric(x<=0)
  else {
    # get the c mixture for the positive lambdas
    beta <- seq(.ruben(lams.pos), min(lams.pos), length = 10)
    c.pos <- NULL
    choice <- 0
    while (is.null(c.pos)) {
      choice <- choice + 1
      c.pos <- .getC(lams.pos, beta[choice], eps)
    }
    beta.pos <- beta[choice]
  
    # get the c mixture for the negative lambdas
    beta <- seq(.ruben(lams.neg), min(lams.neg), length = 10)
    c.neg <- NULL
    choice <- 0
    while (is.null(c.neg)) {
      choice <- choice + 1
      c.neg <- .getC(lams.neg, beta[choice], eps)
    }
    beta.neg <- beta[choice]
                 
    # find the p-value
    i <- 1:length(c.pos)
    df.1 <- m.pos+2*(i-1)
    pfs <- sapply(1:length(c.neg), function(j) {
      df.2 <- m.neg+2*(j-1)
      qq <- (x*beta.neg/beta.pos)*(df.2/df.1)
      pf(q=qq, df1=df.1, df2=df.2, lower.tail=FALSE)
    })

    p.value <- drop(crossprod(c.pos, pfs) %*% c.neg) 
  }

  return(p.value)           
}


####################################
# The Imhof method for x=0 and central variables only.
####################################
.pImhof <- function(lams, eps = 1e-10) {
  lams <- lams[lams != 0]
  integrand <- function(u) {          # the Imhof integrand. Domain: 0...Inf
    u <- u
    theta <- 0.5 * colSums(atan(outer(lams,u))) # - 0.5 * x * u
    rho <- exp(colSums(0.25 * log(1 + outer(lams^2,u^2))))
    out <- ifelse(u==0, sum(lams)/2, sin(theta)/(u*rho))
    out
  }
  tr.integrand <- function(v) {       # Transformation of the integrand. Domain: 0...1
    K <- sum(abs(lams))/20            # Scaling constant of the transformation (to make it invariant to the scale of lams)
    0.5 + integrand(-log(1-v)/K) / (pi*K*(1-v))
  }
  rt <- max(50*.Machine$double.eps, eps)
  res <- try(integrate(tr.integrand, 0, 1, rel.tol = rt), silent=TRUE)
  if (is(res, "try-error")) {
    out <- list(value = NA, error = NA)
  } else {
    out <- list(value = res$value, error = res$abs)
  }
  out
}
  