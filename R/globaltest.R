require("methods")

#==========================================================
# Function "globaltest" performs the global test on (a list of) 
#    subsets of the data
#
# X is a data matrix (n rows are samples, p columns are genes) 
#    or Biobase exprSet.
# NB: if the dimension of X does not fit the dimension of Y, 
#    but t(X) does, t(X) is used.
# Code missing values as NA
#
# Y is a vector of n clinical outcomes or the name or index of one of 
#    the phenoData variables ( if X is an exprSet) 
#
# test.genes.lst = list of vectors
#       vector can be a length p vector of 0 and 1
#           (1 = selected; 0 = not selected).
#       or a vector with the numbers of the selected genes
#           e.g. c(5,8,9) to test genes 5, 8 and 9.
#       or a vector with the ids in the rownames of X
#           e.g c("AA173143","AA461092","AA487442","AA490243","R26706","R61845")
# 
# OPTIONS:
# model = "logistic": (default).
# model = "linear": for continuous Y.
# levels: vector of groups to test in Y. Not needed if Y is binominal
#           If Y contains > 2 levels then the following methods are used
#           if levels contains 1 value: this groups is tested against all other samples
#           if levels contains 2 values: these groups are tested against each other
#           TODO: if levels contain > 2 values: test all groups
# adjust = data.frame: confounders or riskfactors for which the test must be adjusted
#       may be coded as names or indices of variables in the PhenoData slot of X (if exprSet)
# permutation = TRUE: use permutation version of the test
# nperm: number of permutations to be used
# sampling = TRUE: compare to sampled pathways of same size
# ndraws = number: number of "pathways" sampled (default = 10^3).
# verbose = TRUE: print progress data
#
# RESULT
# array with 7 columns containing
# p.value, Q, EQ, seQ, comparative.p, length of vector, rows found in X
# where comparative.p is the fraction of random pathways of the same size 
#   with a lower p-value
# TODO: if #levels > 2 then the result is a list of arrays
#==========================================================

globaltest <- function(X, Y, test.genes = NULL,
                        model = NULL,
                        levels = NULL,
                        adjust = NULL,
                        permutation = FALSE,
                        nperm = NULL,
                        sampling = FALSE,
                        ndraws = NULL,
                        verbose = TRUE)

{
    # check for correct input of X and Y:
    # 1: coerce Y into a vector
    if ( is(X, "exprSet") & length(Y)==1 ) {
      if ((Y %in% names(pData(X))) | (Y %in% 1:ncol(pData(X)))) {
        Y <- pData(X)[,Y]
        names(Y) <- sampleNames(X)
      } else { 
        stop("Y not among pData(X) variables", call. = FALSE)
      }
    } else {
      if ( !is.vector(Y) & !is.matrix(Y) & !is.data.frame(Y) )
          stop("Y should be of type 'vector'", call. = FALSE)
      if ( is.data.frame(Y) )
          Y <- as.matrix(Y)
      if ( is.matrix(Y) ) {
        if ( nrow(Y) == 1 )
          Y <- t(Y)
        if ( ncol(Y) == 1 ) {
          samplenamesY <- rownames(Y)
          Y <- as.vector(Y)
          names(Y) <- samplenamesY
        } else {
          stop("Y should be a single vector.", call. = FALSE)
        }
      }
    }
    
    # 2: Determine the model from the input Y
    if (is.null(model)) {
      if (is.factor(Y) | (length(levels(factor(Y))) <= 2) | !is.null(levels)) {
        model <- 'logistic'
      } else {
        if (is.numeric(Y)) {
          model <- 'linear'
        } else {
          stop("model could not be determined from input Y: please specify model", call. = FALSE)
        }
      } 
    }
    model<- match.arg(model, c('linear','logistic')) 

    # 3: coerce X into a matrix
    if (is(X, "exprSet")) {
        exX <- X
        X <- exprs(X)
    } else {
        exX <- NULL
    }
    if (!is.vector(X) & !is.data.frame(X) & !is.matrix(X))
        stop("X should be of type 'matrix'.", call. = FALSE)
    X <- as.matrix(X)

    # 4: Preparation of X and Y in the logistic model using option levels
    if (model=='logistic') {
      if (is.null(levels)){
        # Only 2 levels should be here, test for 1 now, later checks will find other errors
        levels <- levels(factor(Y))
        if (length(levels)==1) 
          stop("There should be more than 1 group in the data.", call. = FALSE)
      } else {
        if (!all(levels %in% levels(factor(Y))))
          stop("input 'levels' incompatible with levels of Y", call. = FALSE)
        }
      if (length(levels)==2 ) {
      # create a subset of samples
        X <- as.matrix(X[,Y == levels[1] | Y == levels[2]])
        exX <- exX[,Y == levels[1] | Y == levels[2]]
        Y <- Y[Y == levels[1] | Y == levels[2]]
      }
      if (length(levels)<=2 ) {
        samplenamesY <- names(Y)
        Y <- (Y != levels[1])
        names(Y) <- samplenamesY
      }else{ 
        stop("No more than 2 groups can be tested. Use option: levels", call. = FALSE)    
      }
      Y <- 0 + Y
    } else {
      levels <- numeric(0)
    }

    # 5: check if X and Y are numeric
    if (!is.numeric(X) | !is.numeric(Y))
        stop("X and Y must be numeric values.", call. = FALSE)
    if (any(is.na(Y)))
        stop("missing values not allowed in Y.", call. = FALSE)
    
    # 6: check dimensions of X and Y 
    n <- length(Y)
    if (ncol(X) == n) {
        # probably all exprsets will be transposed here
        p <- nrow(X)
        X <- t(X)
    }else{
        if (nrow(X) == n)
             p <- ncol(X)
        else
            stop("Dimensions of X and Y don't match.", call. = FALSE)
    }
    if (n == p)
        warning("As many samples as genes.\n Columns are assumed samples; rows assumed genes.", call.=FALSE)

    # 7: extract genenames and compare samplenames
    genenames <- colnames(X)
    if (  (is.null(rownames(X))) & (!is.null(names(Y))) ) {
      rownames(X) <- names(Y)
    }else{
      if (  (!is.null(rownames(X))) & (is.null( names(Y) )) ) {
        names(Y) <- rownames(X)
      }else{ 
        if ( !all(names(Y) == rownames(X)) )
          stop("Sample names in X inconsistent with sample names in Y.", call. = FALSE)
      }
    }

    # Check for correct input of adjust
    if (!is.null(adjust)) { 
      if (permutation | !is.null(nperm))
        stop("permutation version not compatible with adjusted test", call. = FALSE)
      ready <- FALSE
      if (!is.null(exX)) {
        if (all(adjust %in% names(pData(exX))) | all(adjust %in% 1:ncol(pData(exX)))) {
          adjustnames <- adjust
          adjustmodel <- as.formula(paste("Y ~ ", paste(adjust, collapse = "+")))
          adjust <- pData(exX)
          adjustsamplenames <- TRUE
          ready <- TRUE
        } else {
          if (any(adjust %in% names(pData(exX))))
            stop("not all variables in adjust found among names(pData(X))", call. = FALSE)
        }
      }
      if (!ready) { 
        if (is.vector(adjust)) {
          adjustsamplenames <- !is.null(names(adjust))
          adjust <- data.frame(adjust)
          adjustnames <- names(adjust)
          adjustmodel <- as.formula(paste("Y ~ ", paste(names(adjust), collapse = "+")))
        } else {
          if (is.matrix(adjust)) {
            adjustsamplenames <- !is.null(rownames(adjust))
            adjust <- data.frame(adjust)
            adjustnames <- names(adjust)
            adjustmodel <- as.formula(paste("Y ~ ", paste(names(adjust), collapse = "+")))
          } else {
            if (is.data.frame(adjust)) {
              adjustsamplenames <- TRUE
              adjustnames <- names(adjust)
              adjustmodel <- as.formula(paste("Y ~ ", paste(names(adjust), collapse = "+")))
            } else {
              stop("adjust should be 'vector', 'matrix' or 'data.frame'", call. = FALSE)
            }
          }
        }
        if (nrow(adjust) != n) 
          stop("the number of rows of adjust does not match the length of Y", call. = FALSE)
      }

      # consistency and missing value check
      if (adjustsamplenames & any(rownames(adjust) != names(Y)))
        stop("samplenames in adjust do not match samplenames in Y or X", call. = FALSE)
      if (any(is.na(adjust[,adjustnames])))  
        stop("missing values not allowed in adjust", call. = FALSE)
      # make the adjustment X-matrix by expanding the factors into dummies
      Z <- matrix(1, n, 1)
      for (name in adjustnames) {
        if (is.numeric(adjust[,name])) {
          Z <- cbind(Z, adjust[,name])
        } else {
          if (is.factor(adjust[,name])) {
            lvls <- levels(adjust[,name])
            dummies <- outer(adjust[,name], lvls, "==")
            levelpresent <- apply(dummies, 2, any)
            dummies <- as.matrix(dummies[,levelpresent])
            if (ncol(dummies) > 1)
              Z <- cbind(Z, dummies[,2:ncol(dummies)])
          } else {
            stop("adjust should contain numeric or categorical data", call. = FALSE)
          }
        }
      }
    } else {
      adjust <- data.frame(Y)
      Z <- matrix(rep(1, times = n), n, 1)
      adjustmodel <- Y ~ 1
    }
    nadjust <- ncol(Z)
    
    # Fit adjustmodel and adjust Y and X
    # Standardize Y in case of the linear model or unadjusted logistic (same variance for i=1,...,n)
    if (model == 'linear') {
      fit <- lm(adjustmodel, adjust)
      old.sumYY <- sum(Y * Y)
      Y <- residuals(fit)
      if (df.residual(fit) > 0)
        mu2 <- sum(Y*Y)/df.residual(fit)
      else
        stop("no degrees of freedom left after adjustment", call. = FALSE)
      Rsquare <- sum(Y * Y) / old.sumYY
      Y <- Y / sqrt(mu2)
      IminH <- diag(n) - Z %*% solve(t(Z) %*% Z, t(Z))
      adjustmatrix <- t(IminH)
    } else {
      if (model == 'logistic') {
        if (nadjust > 1) {
          old.mu <- mean(Y)
          fit <- glm(adjustmodel, adjust, family = binomial)
          mu <- fitted.values(fit)
          Y <- Y - mu
          mu2 <- mu * (1 - mu)
          Rsquare <- sum(mu2) / (n * old.mu * (1 - old.mu))
          IminH <- diag(n) - diag(mu2) %*% Z %*% solve(t(Z) %*% diag(mu2) %*% Z, t(Z))
          adjustmatrix <- diag(sqrt(mu2)) %*% t(IminH)
        } else {
          mu <- rep(mean(Y), times = n)
          mu2 <- mu * (1 - mu)
          Y <- (Y - mu) / sqrt(mu2)
          Rsquare <- 1
          IminH <- diag(n) - matrix(rep(1/n, times = n*n), n, n)
          adjustmatrix <- t(IminH)
        }
      }
    }
    if (Rsquare < 0.01)
      stop("not enough variance in Y remaining after adjustment", call. = FALSE)
    
    # Impute missing values of X and rescale X
    col.meanX <- colMeans(X, na.rm = TRUE)
    X <- X - rep(1, times=n) %o% col.meanX
    X[is.na(X)] <- 0
    adjX <- adjustmatrix %*% X
    Einf <- sqrt(sum(adjX * adjX) / (10*p))
    X <- X / Einf
    
    # check for correct input of sampling and ndraws
    if ( (!is.logical(sampling)) | (length(sampling) != 1) )
      stop("Option sampling should be either TRUE or FALSE", call. = FALSE)
    if (!is.null(ndraws)){
      if ( !(ndraws == as.integer(ndraws)) | (ndraws < 1) | (length(ndraws) != 1) ) {
        stop("Option ndraws should be a single positive integer", call. = FALSE)
      }else{
         sampling <- TRUE
      }
    }else{
      if (sampling)
        ndraws <- 10^3
    }        
    
    # check for correct input of permutation and nperm
    if ( (!is.logical(permutation)) | (length(permutation) != 1) )
        stop("Option permutation should be either TRUE or FALSE", call. = FALSE)
    if (!is.null(nperm)){
        if ( !(nperm == as.integer(nperm)) | (nperm < 1) | (length(nperm) != 1) ) {
            stop("Option nperm should be a single positive integer", call. = FALSE)
        }else{
            permutation <- TRUE
        }
    }else{
        if (permutation)
            nperm <- 10^4
    }        
 
    # check incompatible demands
    if (sampling & permutation) 
        stop("sampling = TRUE and permutation = TRUE not allowed at the same time")
   
    # coerce test.genes into a list
    if ( is.list(test.genes) ) {
        test.genes.lst <- test.genes
    }else{
        test.genes.lst <- list(test.genes)
    }

    # prepare output per pathway
    res <- matrix(0, nrow = length(test.genes.lst), ncol = 7)
    colnames(res) <- c("path.n","test.n","Q","EQ","seQ","p.val","comp.p")
    rownames(res) <- names(test.genes.lst)

    verbose <- ( verbose & (sampling | permutation) & (length(test.genes.lst) > 20) )
    if (verbose) {
      cat(length(test.genes.lst), "genesets to be processed. Each dot represents", 
        length(test.genes.lst) %/% 20, "pathways.\n")
      countprogress <- 0
    }

    sampledist<-list(NULL)
    
    # Repeat the following for each geneset
    for (index in 1:length(test.genes.lst)) {
        
      if (verbose) {
        countprogress <- countprogress + 1
        if ( countprogress %% (length(test.genes.lst) %/% 20) == 0) {
          cat(".")
        }
      }
      test.genes <- test.genes.lst[[index]] 
      # check for correct input of test.genes
      # 1: set default
      if (is.null(test.genes)) {
        if (!is.null(genenames))
          test.genes <- genenames
        else
          test.genes <- 1:p
        res[index, "path.n"] <- p
      } else {
        # 2: coerce test.genes into a vector
        if ( !is.vector(test.genes) )
          stop("test.genes should be a (list of) vector(s)", call. = FALSE)
      }
  
      # 3: check format of test.genes
      if ( all(test.genes %in% c(0,1)) & (length(test.genes) == p) ) {
        test.genes <- (test.genes == 1)
        test.names <- names(test.genes)[test.genes]
        res[index, "path.n"] <- sum(test.genes)
      } else {
        if (is.character(test.genes))  {
          res[index, "path.n"] <- length(test.genes)
          test.genes <- intersect(test.genes, genenames)
          test.names <- NULL
        } else {  
          if ( all(test.genes %in% 1:p) ) {
            test.names <- names(test.genes)
            res[index, "path.n"] <- length(test.genes)
          } else {
            stop("option 'test.genes' should be a (list of) vector(s) of gene names or numbers.", call. = FALSE)
          }
        } 
      }  
        
      # 4: check compatibility of names
      if (!is.character(test.genes)) {
        if (  !is.null(genenames )  &  !is.null(test.names) ) {
          if( any( genenames[test.genes] != test.names ) ) 
            warning("Gene names in X inconsistent with gene names in test.genes.", call. = FALSE)
        }
      }
    
      # select genes to be tested
      X.sel <- as.matrix(X[,test.genes])
      # Number of selected genes that are in array (test.genes can be larger than that (in res[index, "path.n"] earlier)
      m <- ncol(X.sel)
      res[index, "test.n"] <- m
      if (m>0) {

        # calculate test statistic Q
        R <- ( t(IminH) %*% X.sel %*% t(X.sel) %*% (IminH) ) / m
        XX <- (X.sel %*% t(X.sel)) / m
        Q <- ( Y %*% XX %*% Y )
         
        # Expectation and variance of Q and p-value
        if (!permutation) {
            # Asymptotic p.value
            EQ <- sum(diag(R))
            trRR <- sum(R*R)
            tr2R <- EQ^2
            trR2 <- sum(diag(R*R))
            if (model == 'logistic') {
              if (nadjust == 1) {
                mu1 <- mu[1] 
                mumu <- mu2[1]
                K <- ( 1 - 6 * mu1 + 6 * mu1^2 ) / mumu
                varQ <- K * ( trR2 - tr2R / n ) + 2 * trRR - 2 * tr2R / (n-1)
              } else {
                # in this case Y is not standardized
                RV <- R %*% diag(mu2)
                EQ <- as.numeric(sum(diag(RV)))
                mu4 <- mu * (1-mu)^4 + (1-mu) * mu^4
                varQ <-  sum(diag(R) * diag(R) * (mu4 - 3 * mu2 * mu2)) + 2 * sum(diag(RV %*% RV))
              }
            }else{ 
              # model = 'linear'
              varQ <- (2 / (n - nadjust + 2)) * ( (n - nadjust) * trRR - tr2R )
            }    
            seQ <- sqrt(varQ)
            scl <- varQ / (2 * EQ)
            dfr <- EQ / scl
            p.value <- pf ( (scl * dfr / Q), 10^10, dfr )
        }else{
          # Permutation p.value
          permQ <- function(stR, YY, npm) {
            pms <- apply( matrix(rnorm(n * npm), n, npm), 2, sort.list )
            YY.pm <-  matrix( YY[pms], n, npm )
            colSums(( stR %*% YY.pm ) * YY.pm)
          }
          Qs <- numeric(nperm)
          chunk <- 5000
          nchunks <- trunc(nperm / chunk)
          rest <- nperm - nchunks * chunk
          anonymousY <- Y
          names(anonymousY) <- NULL
          if (nchunks > 0)
            for (i in 1:nchunks) {
              Qs[(chunk * (i-1) + 1):(chunk * i)] <- permQ( XX, anonymousY, chunk ) 
              NULL
            }
          if (rest > 0)
            Qs[(chunk * nchunks + 1):nperm] <- permQ( XX, anonymousY, rest ) 
          p.value <- sum(Qs >= rep(Q, times=nperm)) / nperm
          EQ <- mean(Qs)
          seQ <- sd(Qs)
        }
 
        if (sampling) {
          # Check if this geneset size has been sampled already
          if (m > length(sampledist)) length(sampledist)<-m
          if (is.null(sampledist[[m]])) {
            ps <- vector("numeric",ndraws)
            for (sample in 1:ndraws){
              X.sample <- as.matrix(X[,sample(1:p,m)])
              R <- ( t(IminH) %*% X.sample %*% t(X.sample) %*% (IminH) ) / m
              XX <- (X.sample %*% t(X.sample)) / m
              Q.sample <- as.numeric(Y %*% XX %*% Y)
              EQ.sample <- sum(diag(R))
              trRR <- sum(R*R)
              tr2R <- EQ.sample*EQ.sample
              trR2 <- sum(diag(R*R))
              if (model == 'logistic') {
                if (nadjust == 1) {
                  mu1 <- mu[1] 
                  mumu <- mu2[1]
                  K <- ( 1 - 6 * mu1 + 6 * mu1^2 ) / mumu
                  varQ <- K * ( trR2 - tr2R / n ) + 2 * trRR - 2 * tr2R / (n-1)
                } else {
                  RV <- R %*% diag(mu2)
                  EQ.sample <- as.numeric(sum(diag(RV)))
                  musq <- mu*mu
                  mu4 <- mu - 4*musq + 6*musq*mu - 3*musq*musq  
                  VarQ <-  sum(diag(R) * diag(R) * (mu4 - 3 * mu2 * mu2)) + 2 * sum(diag(RV %*% RV))
                }
              }else{ 
                # model = 'linear'
                varQ <- (2 / (n - nadjust + 2)) * ( (n - nadjust) * trRR - tr2R )
            }    
            scl <- varQ / (2 * EQ.sample)
            dfr <- EQ.sample / scl
            ps[sample] <- pf ( (scl * dfr / Q.sample), 10^10, dfr )
            }
            sampledist[[m]]<-ps
          }
          comparative.p <- mean(sampledist[[m]] < (rep(p.value, times=ndraws)) * (1+10^-6))
        }else 
          comparative.p <- NA
        
        # the returns for this geneset
        test.genes.lst[[index]] <- test.genes
        res[index, "p.val"] <- p.value
        res[index, "Q"] <- Q
        res[index, "EQ"] <- EQ
        res[index, "seQ"] <- seQ
        res[index, "comp.p"] <- comparative.p
      }else{
        # returns for an empty geneset
        res[index, "p.val"] <- NA
        res[index, "Q"] <- NA
        res[index, "EQ"] <- NA
        res[index, "seQ"] <- NA
        res[index, "comp.p"] <- NA
      }
  }
  
  # returns an object of type "gt.result"
  new("gt.result",
    res = res, 
    X = X,
    Y = Y,
    test.genes = test.genes.lst,
    adjustmatrix = adjustmatrix,
    Rsquare = Rsquare,
    model = model,
    levels = levels,
    df.adjust = nadjust)
}
#==========================================================


#==========================================================
setClass("gt.result", representation(
    res = "matrix", 
    X = "matrix",
    Y = "vector",
    test.genes = "list",
    adjustmatrix = "matrix",
    Rsquare = "numeric",
    model = "character",
    levels = "vector",
    df.adjust = "integer")
)

#==========================================================
# Function "show" prints out a result of type "gt.result"
# such as results from a call to "globaltest"
#==========================================================

setMethod("show", "gt.result", function(object)
{
  npathways <- length(object@test.genes)
  nsamples <- length(object@Y)
  ngenes <- dim(object@X)[2]

  cat("Global Test result:\n")
  if (npathways == 1)
    cat("Data:", nsamples, "samples with", ngenes, "genes; 1 pathway tested\n")
  else
    cat("Data:", nsamples, "samples with", ngenes, "genes;", npathways, "pathways tested\n")
  cat("Model:", object@model, "\n")
  if (object@Rsquare != 1)
    cat("Adjusted:", 100 * round(object@Rsquare, 3), "% of variance of Y remains after adjustment\n")
  cat("\n")
    
  res <- data.frame(object@res)
  if ( all(is.na(res[,"comp.p"])) ) {
    res <- res[,1:6]
    colnames(res) <- c("genes","tested","Statistic Q","Expected Q","sd of Q","p-value")
  }else{
    colnames(res) <- c("genes","tested","Statistic Q","Expected Q","sd of Q","p-value","comp. p")
  }
  print(signif(res, digits = 5))
})


#==========================================================
# Two functions to extract relevant information from 
# a gt.result object
#==========================================================
if( !isGeneric("result") )
    setGeneric("result", function(object) standardGeneric("result"))

setMethod("result", "gt.result",
            function(object) 
{
    res <- object@res
 
    if ( all(is.na(res[,"comp.p"])) ) {
      res <- res[,c("path.n","test.n","Q","EQ","seQ","p.val")]
    }else{
      res <- res[,c("path.n","test.n","Q","EQ","seQ","p.val","comp.p")]
    }
    
    res
})

#==========================================================
if( !isGeneric("p.value") )
    setGeneric("p.value", function(gt) standardGeneric("p.value"))

setMethod("p.value", "gt.result",
            function(gt) gt@res[,"p.val"])


#==========================================================
# Function "permutations" compares the theoretical values 
# of EQ, seQ and the p.value to values based on permutations 
# of the clinical outcome, which may be better for small 
# sample sizes. It summarizes the Q-values for the permutations 
# in a histogram, in which an arrow points out the value of the 
# true Q, for comparison
#==========================================================

permutations <- function(gt, geneset = NULL, nperm = 10^4)
{
    # check correct input of gt
    if ( !is(gt, "gt.result"))
      stop("permutations should be applied to a globaltest result", call. = FALSE)
    if (gt@df.adjust > 1)
      stop("the permutation procedure is not applicable for the adjusted global test", call. = FALSE)      
        
    # check correct input of nperm
    if ( !(nperm > 0) )
      stop("option nperm should be a positive integer", call. = FALSE)
        
    # extract the right test.genes vector
    test.genes.lst <- gt@test.genes
    if ( is.null(geneset) ) {
        if ( (length(test.genes.lst) == 1) ) {
            test.genes <- test.genes.lst[[1]]
            geneset <- 1
        }else{
            stop("option geneset is missing with no default")
        }
    }else{
        if ( is.character(geneset) ) {
            if ( length(intersect(names(test.genes.lst), geneset )) == 1 ) {
                test.genes <- test.genes.lst[[geneset]]
            }else{
                ifelse( (length(intersect(names(test.genes.lst), geneset )) == 0),
                    stop("requested geneset was not among the tested ones", call. = FALSE),
                    stop("more than one geneset given", call. = FALSE) )
            }
        }else{
            if ( length(intersect(1:length(test.genes.lst), geneset)) == 1 ) {
                test.genes <- test.genes.lst[[geneset]]
            }else{
                stop("incorrect input of geneset", call. = FALSE)
            }
        }
    }
    if (is.character(test.genes))
      test.genes <- intersect(test.genes, colnames(gt@X))
    if (gt@res[geneset, "test.n"] == 0)
      stop("empty pathway", call. = FALSE)
            
    # Recreate Y and XX
    X <- as.matrix(gt@X[,test.genes])
    m <- ncol(X)
    XX <- (X %*% t(X)) / m
    Y <- gt@Y
    n <- length(Y)
    Q <- gt@res[geneset, "Q"]

    # Recalculate the Q-value for permutations of Y
    permQ <- function(stR, YY, nn, npm) {
      pms <- apply( matrix(rnorm(nn * npm), nn, npm), 2, sort.list )
      YY.pm <-  matrix( YY[pms], nn, npm )
      colSums(( stR %*% YY.pm ) * YY.pm)
    }

    Qs <- numeric(nperm)
    chunk <- 5000
    nchunks <- trunc(nperm / chunk)
    rest <- nperm - nchunks * chunk
    if (nchunks > 0)
      for (i in 1:nchunks){
        Qs[(chunk * (i-1) + 1):(chunk * i)] <- permQ( XX, Y, n, chunk ) 
        NULL
      }
    if (rest > 0)
      Qs[(chunk * nchunks + 1):nperm] <- permQ( XX, Y, n, rest ) 

    # Draw histogram
    hst <- hist(Qs, xlim = c(0, 1.1 * max( c( Qs, Q ) ) ), 
      main = paste( "Histogram of Q for", nperm, "permutations of Y" ),
      xlab = "Values of Q for permuted Y")
    h <- max(hst$counts)
    arrows( Q, h/5, Q, 0 )
    text( Q, h/4.5, 'Q' )

    # No output
    invisible (NULL)
}
#==========================================================


#==========================================================
# Geneplot plots the influence of each gene on the outcome 
#   of the test statistic
# See help(geneplot) for details
#==========================================================

geneplot <- function(gt, geneset = NULL, genesubset = NULL, drawlabels = TRUE, labelsize = 0.6, ...)
{
    # check correct input of gt
    if ( !is(gt, "gt.result"))
        stop("geneplot should be applied to a globaltest result", call. = FALSE)
        
    # extract the right test.genes vector
    test.genes.lst <- gt@test.genes
    if ( is.null(geneset) ) {
        if ( (length(test.genes.lst) == 1) ) {
            test.genes <- test.genes.lst[[1]]
            geneset <- 1
        }else{
            stop("option geneset is missing with no default")
        }
    }else{
        if ( is.character(geneset) ) {
            if ( length(intersect(names(test.genes.lst), geneset )) == 1 ) {
              test.genes <- test.genes.lst[[geneset]]
            }else{
              ifelse( (length(intersect(names(test.genes.lst), geneset )) == 0),
                stop("requested geneset was not among the tested ones", call. = FALSE),
                stop("more than one geneset given", call. = FALSE) )
            }
        }else{
            if ( length(intersect(1:length(test.genes.lst), geneset)) == 1 ) {
                test.genes <- test.genes.lst[[geneset]]
            }else{
                stop("incorrect input of geneset", call. = FALSE)
            }
        }
    }
    if ( !is.null(genesubset) ) {
        if (all(genesubset %in% test.genes)) {
            test.genes <- intersect(genesubset, test.genes)
        }else{
          if (all(genesubset %in% 1:length(test.genes)))
            test.genes <- test.genes[genesubset]
          else
            stop("genesubset is not a subset of the selected geneset")
        }
    }
    if (is.character(test.genes))
      test.genes <- intersect(test.genes, colnames(gt@X))
    if (gt@res[geneset, "test.n"] == 0)
      stop("empty pathway", call. = FALSE)

    # calculate influence per gene and expected influence
    X <- as.matrix(gt@X[,test.genes])
    m <- ncol(X)
    Y <- gt@Y
    n <- length(Y)
    XY <- t(X) %*% Y  
    influence <- XY * XY
    up <- (sign(XY) == 1)
    adjX <- gt@adjustmatrix %*% X
    Einf <- colSums(adjX * adjX)
    trRR <- Einf * Einf
    tr2R <- Einf * Einf
    trR2 <- colSums(adjX * adjX * adjX * adjX)
    if (gt@model == 'logistic') {
      if (gt@df.adjust == 1) {
        mu1 <- mean(Y < 0) 
        mumu <- mu1 * (1 - mu1)
        K <- ( 1 - 6 * mu1 + 6 * mu1^2 ) / mumu
        varinf <- K * ( trR2 - tr2R / n ) + 2 * trRR - 2 * tr2R / (n-1)
      } else {
        # in this case Y is not standardized
        mu <- abs(Y)
        mu2 <- mu * (1 - mu)
        mu4 <- mu * (1-mu)^4 + (1-mu) * mu^4
        varinf <- t(adjX * adjX * adjX * adjX) %*% (mu4/(mu2 * mu2) - 3) + 2 * trRR
      }
    }else{ 
      # model = 'linear'
      varinf <- (2 / (n - gt@df.adjust + 2)) * ( (n - gt@df.adjust) * trRR - tr2R )
    }   
    sd.inf <- sqrt(varinf)
    # Output: reference vector with gene labels if available
    legend <- colnames(X)
    if (!is.null(legend))
      names(legend) <- 1:length(legend)
    
    if (is.null(labelsize))
      labelsize<-par("cex.axis")
    # Plot reference line and influence per gene and color for up/down regulation
    if (drawlabels & !is.null(legend)){
    # check for space in margin of plot
      plot.new()  
      labwidth<-max(strwidth(legend,"inches",labelsize))
      margins<-par("mai")
      par(new=TRUE,"mai"=c(max(margins[1],labwidth*1.3),margins[2:4]))
      plot( 0, xlim = c(1/2, m+1/2), ylim = c(0, 1.2 * max(influence, Einf) ), col = 0, xlab = "", ylab = "influence", xaxt="n",...)
      axis(1, 1:length(legend), legend, cex.axis = labelsize, las=2)
      par("mai"=margins)
    }else
      plot( 0, xlim = c(1/2, m+1/2), ylim = c(0, 1.2 * max(influence, Einf) ), col = 0, xlab = "genenr", ylab = "influence", ...)
    if (m <= 300) {
        rect(xleft = 1:m - 0.4, xright = 1:m + 0.4, ybottom = rep(0,times=m), ytop = influence, col = (up+2), border=0 )
        nlines <- floor((influence - Einf) / sd.inf)
        nlines[(nlines < 0) | (sd.inf < 10^-3)] <- 0
        for (index in 1:max(nlines)) {
          high <- Einf + apply(cbind(index, nlines),1, min) * sd.inf
          segments( 1:m - 0.4, high, 1:m + 0.4, high )
        }
    } else {
        lines(1:m, influence, lwd = 600 / m, type = 'h', col = (up+2))
    }
    lines( (1:(m+1))-1/2, c(Einf, Einf[m]), type = "s" )
    if (gt@model == 'linear')
      if (may.permute)
        legend(1/2, 1.2 * max(c(influence,Einf)),
          c("positive correlation with Y", "negative correlation with Y"), fil = c(3,2))
      else
        legend(1, 1.2 * max(c(influence,Einf)),
          c("positive correlation with residuals", "negative correlation with residuals"), fil = c(3,2))
    if (gt@model == 'logistic')
      legend(1/2, 1.2 * max(c(influence,Einf)),
        c(paste("higher expression in", gt@levels[2], "samples"), 
        paste("higher expression in", gt@levels[1], "samples")), fil = c(3,2))
       
    invisible(legend)
}
#==========================================================


#==========================================================
# Sampleplot plots the influence of each sample on the outcome 
#   of the test statistic
# See help(sampleplot) for details
#==========================================================
sampleplot <- function(gt, geneset = NULL, samplesubset = NULL, drawlabels = TRUE, labelsize = 0.6,...)
{
    # check correct input of gt
    if ( !is(gt, "gt.result"))
        stop("geneplot should be applied to a globaltest result", call. = FALSE)
        
    # extract the right test.genes vector
    test.genes.lst <- gt@test.genes
    if ( is.null(geneset) ) {
        if ( (length(test.genes.lst) == 1) ) {
            test.genes <- test.genes.lst[[1]]
            geneset <- 1
        }else{
            stop("option geneset is missing with no default")
        }
    }else{
        if ( is.character(geneset) ) {
            if ( length(intersect(names(test.genes.lst), geneset )) == 1 ) {
              test.genes <- test.genes.lst[[geneset]]
            }else{
              ifelse( (length(intersect(names(test.genes.lst), geneset )) == 0),
                stop("requested geneset was not among the tested ones", call. = FALSE),
                stop("more than one geneset given", call. = FALSE) )
            }
        }else{
            if ( length(intersect(1:length(test.genes.lst), geneset)) == 1 ) {
                test.genes <- test.genes.lst[[geneset]]
            }else{
                stop("incorrect input of geneset", call. = FALSE)
            }
        }
    }
    if (is.character(test.genes))
      test.genes <- intersect(test.genes, colnames(gt@X))
    if (gt@res[geneset, "test.n"] == 0)
      stop("empty pathway", call. = FALSE)
    if (!is.null(samplesubset) & !is.vector(samplesubset))
      stop("samplesubset should be of type 'vector'", call. = FALSE)

    # calculate influence per sample and expected influence
    X <- as.matrix(gt@X[,test.genes])
    n <- nrow(X)
    if (!is.null(samplesubset) & !all(samplesubset %in% 1:n) & !all(samplesubset %in% rownames(X)))
      stop("samplesubset should contain names or numbers of samples", call. = FALSE)
    if (is.null(samplesubset))
      samplesubset <- 1:n
    else 
      if (all(samplesubset %in% rownames(X)))
        samplesubset <- match(samplesubset, rownames(X))
    m <- ncol(X)
    Y <- gt@Y
    XXY <- X %*% t(X) %*% Y
    influence <- (n/m) * Y * XXY
    up <- (sign(Y) == 1) 
    if ((gt@model == 'logistic') & (gt@df.adjust > 1)) {
      mu2 <- abs(Y) * (1 - abs(Y))
      tIminH <- diag(1/sqrt(mu2)) %*% gt@adjustmatrix #adjustmatrix <- diag(sqrt(mu2)) %*% t(IminH)
      adjX <- tIminH %*% X
      R <- adjX %*% t(adjX)
      RR <- (R * R) %*% diag(mu2)
    } else {
      tIminH <- gt@adjustmatrix
      adjX <- tIminH %*% X
      R <- adjX %*% t(adjX)
      RR <- R * R
    }
    Einf <- (n/m) * Y * Y * rowSums(adjX*adjX)
    RR <- matrix(t(RR)[diag(n) == 0], n, n-1, byrow = TRUE)
    varinf <- (n/m)^2 * Y * Y * rowSums(RR)
    sd.inf <- sqrt(varinf)
    influence <- influence[samplesubset]
    Einf <- Einf[samplesubset]
    sd.inf <- sd.inf[samplesubset]
    up <- up[samplesubset]
    n <- length(samplesubset)

    # Plot reference line and influence per gene and color for up/down regulation
    rangebars <- max(influence - Einf) - min(influence - Einf)
    minplot <- min(influence - Einf) 
    maxplot <- max(influence - Einf) + 0.2 * rangebars
    # Output: reference vector with sample labels if available
    legend <- rownames(X)[samplesubset]
    if (!is.null(legend))
      names(legend) <- 1:length(legend)
    if (is.null(labelsize))
      labelsize<-par("cex.axis")
    if (drawlabels & !is.null(legend)){
      # check for space in margin of plot
      plot.new()  
      labwidth<-max(strwidth(legend,"inches",labelsize))
      margins<-par("mai")
      par(new=TRUE,"mai"=c(max(margins[1],labwidth*1.3),margins[2:4]))
      plot( 0, xlim = c(1/2, n+1/2), ylim = c(minplot, maxplot), col = 0, xlab = "", ylab = "influence" , xaxt="n",...)
      axis(1, 1:length(legend), legend, cex.axis = labelsize, las=2)
      par("mai"=margins)
    }else
      plot( 0, xlim = c(1/2, n+1/2), ylim = c(minplot, maxplot), col = 0, xlab = "samplenr", ylab = "influence", ...)
    if (n <= 300) {
        rect(xleft = 1:n - 0.4, xright = 1:n + 0.4, ybottom = rep(0,times=n), ytop = influence - Einf, col = (up+2), border=0 )
        nlines <- trunc((influence - Einf) / sd.inf)
        nlines[sd.inf < 10^-3] <- 0
        for (index in 1:max(abs(nlines))) {
          high <- sign(nlines) * apply(cbind(abs(index), abs(nlines)),1, min) * sd.inf
          segments( 1:n - 0.4, high, 1:n + 0.4, high)
        }
    }else{
        lines(1:n, influence - Einf, lwd = 600 / n, type = 'h', col = (up+2))
    }
    segments(1/2, 0, n+1/2, 0)
    if (gt@model == 'linear')
      legend(1/2, maxplot,
        c("positive residual Y", "negative residual Y"), fil = c(3,2))
    if (gt@model == 'logistic')
      legend(1/2, maxplot,
        c(paste(gt@levels[2], "samples"), paste(gt@levels[1], "samples")), fil = c(3,2)) 
        
    invisible(legend)
}
#==========================================================

#==========================================================
# The function regressionplot allows the evaluation of 
#   possibly outlying samples.
# See help(regressionplot) for details
#==========================================================

regressionplot <- function(gt, geneset = NULL, sampleid = NULL,...)
{
    # check correct input of gt
    if (!is(gt, "gt.result"))
        stop("regressionplot should be applied to a globaltest result", call. = FALSE)
    if (gt@df.adjust > 1)
      stop("regressionplot not possible for adjusted globaltest result", call. = FALSE)
        
    # extract the right test.genes vector
    test.genes.lst <- gt@test.genes
    if ( is.null(geneset) ) {
        if ( (length(test.genes.lst) == 1) ) {
            test.genes <- test.genes.lst[[1]]
            geneset <- 1
        }else{
            stop("option geneset is missing with no default")
        }
    }else{
        if ( is.character(geneset) ) {
            if ( length(intersect(names(test.genes.lst), geneset )) == 1 ) {
                test.genes <- test.genes.lst[[geneset]]
            }else{
                ifelse( (length(intersect(names(test.genes.lst), geneset )) == 0),
                    stop("requested geneset was not among the tested ones", call. = FALSE),
                    stop("more than one geneset given", call. = FALSE) )
            }
        }else{
            if ( length(intersect(1:length(test.genes.lst), geneset)) == 1 ) {
                test.genes <- test.genes.lst[[geneset]]
            }else{
                stop("incorrect input of geneset", call. = FALSE)
            }
        }
    }
    if (gt@res[geneset, "test.n"] == 0)
      stop("empty pathway", call. = FALSE)
    if (is.character(test.genes))
      test.genes <- intersect(test.genes, colnames(gt@X))


    # recreate matrix R and find S = Y %o% Y
    X <- as.matrix(gt@X[,test.genes])
    Y <- gt@Y
    R <- as.vector(X %*% t(X))
    S <- as.vector(Y %o% Y)
    n <- length(Y)
    
    # Check correct input of sampleid
    if ( is.null(sampleid) ) {
        sampleid <- 1:n
    }else{
        if ( !xor( (length(intersect(1:n, sampleid)) == 0), (length(intersect(names(Y), sampleid)) == 0) ) )
            stop("Option samplenr incorrect", call. = FALSE)
    }
            
    # Extract relevant entries from S and R
    samples <- rep(FALSE,times = n)
    names(samples) <- rownames(X)
    samples[sampleid] <- TRUE
    samples <- outer(samples, samples, "|")
    lowertriangle <- outer(1:n, 1:n, ">")
    selection <- as.vector(lowertriangle)
    subselection <- as.vector( lowertriangle & samples )
    Rall <- R[selection]
    Sall <- S[selection]
    Rsub <- R[subselection]
    Ssub <- S[subselection]
    Rrest <- R[selection & !subselection]
    Srest <- S[selection & !subselection]
    
    # Draw the plots
    plot(Sall, Rall, xlab = "Covariance between outcomes", ylab = "Covariance between expression profiles", col = 0,...)
    if (length(Rrest) > 0){
        points(Srest, Rrest,col = 4)
        abline(lm(Rall ~ Sall), col = 4)
    }
    points(Ssub, Rsub, col = 2, pch = 4, cex = 1.5)
    abline(lm(Rsub ~ Ssub), col = 2)
    
    # Some explanation
    if ( ( !all(samples) ) & ( !is.null(names(Y)) ) )
        cat("Samples investigated:\n", names(Y)[samplenr], "\n")
         
    # No output
    invisible(NULL)
}   
#==========================================================


#==========================================================
# Function Checkerboard visualizes a result of globaltest
#   by visualizing the covariance matrix between the different samples
# See help(checkerboard) for details
#==========================================================

checkerboard <- function(gt, geneset = NULL, sort = TRUE, drawlabels = TRUE, labelsize = 0.6,...)
{   
    # check correct input of gt
    if ( !is(gt, "gt.result"))
        stop("checkerboard should be applied to a globaltest result", call. = FALSE)
        
    # extract the right test.genes vector
    test.genes.lst <- gt@test.genes
    if ( is.null(geneset) ) {
        if ( (length(test.genes.lst) == 1) ) {
            test.genes <- test.genes.lst[[1]]
            geneset <- 1
        }else{
            stop("option geneset is missing with no default")
        }
    }else{
        if ( is.character(geneset) ) {
            if ( length(intersect(names(test.genes.lst), geneset )) == 1 ) {
                test.genes <- test.genes.lst[[geneset]]
            }else{
                ifelse( (length(intersect(names(test.genes.lst), geneset )) == 0),
                    stop("requested geneset was not among the tested ones", call. = FALSE),
                    stop("more than one geneset given", call. = FALSE) )
            }
        }else{
            if ( length(intersect(1:length(test.genes.lst), geneset)) == 1 ) {
                test.genes <- test.genes.lst[[geneset]]
            }else{
                stop("incorrect input of geneset", call. = FALSE)
            }
        }
    }
    if (gt@res[geneset, "test.n"] == 0)
      stop("empty pathway", call. = FALSE)
    if (is.character(test.genes))
      test.genes <- intersect(test.genes, colnames(gt@X))
    if (!is.logical(sort))
        stop("The option 'sort' should be either TRUE or FALSE", call. = FALSE)

    # recreate matrix R
    X <- as.matrix(gt@X[,test.genes])
    Y <- gt@Y
    R <- X %*% t(X)
        
    # Sort Y if needed and return new samplenrs
    n <- length(Y)
    perm <- sort.list(Y)
    rperm <- sort.list(Y[n:1])
    if ( any(perm != 1:n ) & any(rperm != 1:n) & (sort)){
        newsamplenrs <- matrix( c(1:n, sort.list( perm )), n, 2 ) 
        label <- "sorted samplenr"
        R <- R[perm,perm]
    }else{
        sort = FALSE
        label <- "samplenr"
        newsamplenrs <- matrix( c(1:n, 1:n), n, 2 ) 
    }    
    colnames(newsamplenrs) <- c("samplnr.old", "samplenr.new")
    rownames(newsamplenrs) <- rownames(X)

    # Calculate median non-diagonal element of R
    lowertriangle <- outer( 1:n, 1:n, ">" )
    med <- median(R[lowertriangle])
    # Draw plot
    if (is.null(labelsize))
      labelsize<-par("cex.axis")
    if (drawlabels & !is.null(rownames(newsamplenrs))){
      legend<-rownames(newsamplenrs)[sort(newsamplenrs[,2],index.return=TRUE)$ix]   
      # check for space in margin of plot
      plot.new()  
      labwidth<-max(strwidth(legend,"inches",labelsize))
      margins<-par("mai")
      par(new=TRUE,"mai"=c(max(margins[1],labwidth*1.3),max(margins[2],labwidth*1.3), margins[3:4]))
      image(x = 1:n, y = 1:n, z = R>med, col = rainbow(2, v = c(0,1), s = c(1,0) ), ylab = "", xlab = "", 
          yaxt="n", xaxt = "n", ...)
      axis(2,1:length(legend), legend, cex.axis=labelsize, las=2)
      axis(1,1:length(legend), legend, cex.axis=labelsize, las=2)
      par("mai"=margins)
    }else
      image(x = 1:n, y = 1:n, z = R>med, col = rainbow(2, v = c(0,1), s = c(1,0) ), xlab = label, ylab = label, 
          lab = c(n,n,50/n), ...)
    par(pty = "s")
    invisible(newsamplenrs)
}
#==========================================================

#==========================================================
# .First.lib is called when the package is loaded
# It initializes the object "gt.result"
#   and the its "show" function.
#==========================================================

.First.lib <- function(libname, pkgname, where)
{
    if (missing(where)) {
        where <- match(paste("package:", pkgname, sep=""), search())
        if(is.na(where)) {
            warning(paste("Not a package name: ",pkgname))
            return()
        }
        where <- pos.to.env(where)
    }


    invisible(NULL)
}
#==========================================================
