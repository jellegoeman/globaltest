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
# model = 'logistic': (default).
# model = 'linear': for continuous Y.
# levels = vector of groups to test in Y. Not needed if Y is binominal
#           If Y contains > 2 levels then the following methods are used
#           if levels contains 1 value: this groups is tested against all other samples
#           if levels contains 2 values: these groups are tested against each other
#           TODO: if levels contain > 2 values: test all groups
# sampling = TRUE: compare to sampled pathways of same size
# ndraws = number: number of "pathways" sampled (default = 10^3).
#
# RESULT
# array with 7 columns containing
# p.value, Q, EQ, seQ, comparative.p, length of vector, rows found in X
# where comparative.p is the fraction of random pathways of the same size 
#   with a lower p-value
# TODO: if #levels > 2 then the result is a list of arrays
#==========================================================

globaltest <- function(X, Y, test.genes = NULL,
                        model = 'logistic',
                        levels = NULL,
                        permutation = FALSE,
                        nperm = NULL,
                        sampling = FALSE,
                        ndraws = NULL,
                        verbose = TRUE)

{
    # check for correct input of X and Y:
    # 1: coerce Y into a vector
    if ( is(X, "exprSet") & length(Y)==1 ) {
      Y <- pData(X)[,Y]
      samplenamesY <- sampleNames(X)
      if (model=='logistic') {
        if (is.null(levels)){
          # Only 2 levels should be here, test for 1 now, later checks will find other errors
          levels<-levels(factor(Y))
          if (length(levels)==1) 
            stop("There should be more than 1 group in the data.", call. = FALSE)
        }
        if (length(levels)==2 ) {
          # create a subset of samples
          samplenamesY <- samplenamesY[Y == levels[1] | Y == levels[2]]
          X <- X[,Y == levels[1] | Y==levels[2]]
          Y <- Y[Y == levels[1] | Y==levels[2]]
        }
        if (length(levels)<=2 ) {
          Y <- (Y == levels[1])
        }else{ 
          stop("No more than 2 groups can be tested.", call. = FALSE)    
        }
      }
    }else{
      if ( !is.vector(Y) & !is.matrix(Y) & !is.data.frame(Y) )
          stop("Y should be of type 'vector'.", call. = FALSE)
      if ( is.data.frame(Y) )
          Y <- as.matrix(Y)
      if ( is.matrix(Y) ) {
        if ( dim(Y)[1] == 1 )
          Y <- t(Y)
        if ( dim(Y)[2] == 1 ) {
          namesY <- rownames(Y)
          Y <- as.vector(Y)
          names(Y) <- namesY
        } else {
          stop("Y should be a single vector.", call. = FALSE)
        }
      }
      samplenamesY <- names(Y)
      names(Y) <- NULL
    }

    # 2: coerce X into a matrix
    if (is(X, "exprSet"))
        X <- exprs(X)
    if (!is.vector(X) & !is.data.frame(X) & !is.matrix(X))
        stop("X should be of type 'matrix'.", call. = FALSE)
    X <- as.matrix(X)

    # 3: check if X and Y are numeric
    if (is.logical(Y))    
        Y <- 0 + Y
    if (!is.numeric(X) | !is.numeric(Y))
        stop("X and Y may only have numeric values.", call. = FALSE)
    if (any(is.na(Y)))
        stop("Missing values not allowed in Y.", call. = FALSE)
    
    # 4: check dimensions of X and Y 
    n <- length(Y)
    pn <- dim(X)    
    if (pn[2] == n){
        # probably all exprsets will be transposed here
        p <- pn[1]
        X <- t(X)
    }else{
        if (pn[1] == n)
             p <- pn[2]
        else
            stop("Dimensions of X and Y don't match.", call. = FALSE)
    }
    if (n == p)
        warning("As many samples as genes.\n Columns are assumed to represent samples; rows assumed genes.", call.=FALSE)
        
    # 5: extract genenames and samplenames
    genenames <- colnames(X)
    samplenamesX <- rownames(X)
    if (  ( is.null(samplenamesX) ) & ( !is.null(samplenamesY) ) ) {
      samplenames <- samplenamesY
    }else{
      if (  ( !is.null(samplenamesX) ) & ( is.null( samplenamesY ) ) ) {
        samplenames <- samplenamesX
      }else{ 
        if ( all(samplenamesY == samplenamesX) ) {
          samplenames <- samplenamesY
        }else{
          stop("Sample names in X inconsistent with sample names in Y.", call. = FALSE)
        }
      }
    }
    rownames(X) <- samplenames


    # Check for correct input of model
    if ( (model != 'linear') & (model != 'logistic') )
        stop("Option model should be either 'linear' or 'logistic'.", call. = FALSE)
    if (model == 'logistic'){
        Y1 <- (Y == max(Y))
        Y0 <- (Y == min(Y))
        if ( ( sum(Y0 | Y1) == n ) & ( max(Y) != min(Y) ) )
            Y <- Y1
        else
            stop("For the logistic model Y should take exactly two values.\n For a continuous Y use option: model = 'linear'.", call. = FALSE)
    }
    
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
        if ( !(nperm == as.integer(ndraws)) | (nperm < 1) | (length(nperm) != 1) ) {
            stop("Option nperm should be a single positive integer", call. = FALSE)
        }else{
            permutation <- TRUE
        }
    }else{
        if (permutation)
            nperm <- 10^4
    }        
 
    if (sampling & permutation) 
        stop("sampling = TRUE and permutation = TRUE not allowed at the same time: this takes too long")
   
    # coerce test.genes into a list
    if ( is.list(test.genes) ) {
        test.genes.lst <- test.genes
    }else{
        test.genes.lst <- list(test.genes)
    }

    # center Y
    mu <- mean(Y)
    Y <- Y - mu

    # prepare output per pathway
    res <- matrix(0, nrow = length(test.genes.lst), ncol = 7)
    colnames(res) <- c("path.n","test.n","Q","EQ","seQ","p.val","comp.p")
    rownames(res) <- names(test.genes.lst)

    # Center genes and set missing values to zero
    col.meanX <- colMeans(X, na.rm = TRUE)
    X <- X - rep(1, times=n) %o% col.meanX
    X[is.na(X)] <- 0

    verbose <- ( verbose & (sampling | permutation) & (length(test.genes.lst) > 20) )
    if (verbose) {
      cat(length(test.genes.lst), "genesets to be processed. Each dot represents", 
        length(test.genes.lst) %/% 20, "pathways.\n")
      countprogress <- 0
    }

    sampledist<-list(NULL)
    
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
          test.genes <- rep(TRUE,times = p)
          if (is.null(genenames)) {
            test.genes.lst[[index]] <- test.genes
          }else{
            test.genes.lst[[index]] <- genenames
          }
          res[index, "path.n"] <- p
      }else{
          # 2: coerce test.genes into a vector
          if ( !is.vector(test.genes) & !is.matrix(test.genes) & !is.data.frame(test.genes) )
              stop("test.genes should be of type 'vector'.", call. = FALSE)
          if ( is.data.frame(test.genes) )
              test.genes <- as.matrix(test.genes)
          if ( is.matrix(test.genes) ) {
              if ( dim(test.genes)[1] == 1 )
                  test.genes <- t(test.genes)
              if ( dim(test.genes)[2] == 1 ) {   
                 genes <- rownames(test.genes)     
                 test.genes <- as.vector(test.genes)
                 names(test.genes) <- genes
              } else
                  stop("test.genes should be a single vector.", call. = FALSE)
          }
  
          # 3: check format of test.genes
          if ( all( (test.genes == 0) | (test.genes == 1) ) & (length(test.genes) == p) & (!all(test.genes==0)) ) {
              test.genes <- (test.genes == 1)
              test.names <- names(test.genes)[test.genes]
              res[index, "path.n"] <- sum(test.genes)
          } else {
            if (is.character(test.genes))  {
              res[index, "path.n"] <- length(test.genes)
              test.genes <- intersect(test.genes, genenames)
              test.names <- test.genes
            } else {  
              if ( length(intersect(1:p, test.genes) == length(test.genes)) ) {
                  test.genes <- sort(test.genes)
                  test.names <- names(test.genes)
                  res[index, "path.n"] <- length(test.genes)
              }
              else
                  stop("Option test.genes should have 0 or 1 for each gene\n or be a list of gene names or numbers.", call. = FALSE)
            } 
          }  
        
          # 4: check compatibility of names
          if (!is.character(test.genes)) {
            if (  !is.null(genenames )  &  !is.null( test.names ) ) {
                 if( any( genenames[test.genes] != test.names ) ) 
                     warning("Gene names in X inconsistent with gene names in test.genes.", call. = FALSE)
            }
          }
      }
    
      # select genes to be tested
      X.cent <- as.matrix(X[,test.genes])
      # Number of selected genes that are in array (test.genes can be larger than that (in res[index, "path.n"] earlier)
      m <- dim(X.cent)[2]
      res[index, "test.n"] <- m
      if (m>0) {

        # calculate test statistic Q
        R <- ( X.cent %*% t(X.cent) ) / m
        mu2 <- switch(model,
            'linear' = var(Y),
            'logistic' = mu * (1-mu))
        Q <- as.numeric( ( Y %*% R %*% Y ) / mu2 )
 
        # Expectation and variance of Q and p-value
        if (!permutation) {
            # Asymptotic p.value
            EQ <- sum(diag(R))
            trRR <- sum(R^2)
            tr2R <- EQ^2
            trR2 <- sum(diag(R^2))
            if (model == 'logistic'){
                K <- ( 1 - 6 * mu + 6 * mu^2 ) / mu2
                varQ <- K * ( trR2 - tr2R / n ) + 2 * trRR - 2 * tr2R / (n-1)
            }else{
                varQ <- (2 / (n+1)) * ( (n-1) * trRR - tr2R )
            }    
            seQ <- sqrt(varQ)
            scl <- varQ / (2 * EQ)
            dfr <- EQ / scl
            p.value <- pf ( (scl * dfr / Q), 10^6, dfr )
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
            for (i in 1:nchunks){
                Qs[(chunk * (i-1) + 1):(chunk * i)] <- permQ( R/mu2, Y, chunk ) 
                NULL
            }
            if (rest > 0)
                Qs[(chunk * nchunks + 1):nperm] <- permQ( R/mu2, Y, rest ) 
            p.value <- sum(Qs >= rep(Q, times=nperm)) / nperm
            EQ <- mean(Qs)
            seQ <- sd(Qs)
        }
 
        if (sampling) {
            # Check if this geneset size has been sampled already
            if (m>length(sampledist)) length(sampledist)<-m
            if (is.null(sampledist[[m]])) {
                ps<- vector("numeric",ndraws)
                for (sample in 1:ndraws){
                    X.sample <- as.matrix(X[,sample(1:p,m)])
                    R.sample <- (X.sample %*% t(X.sample)) / m
                    Q.sample <- as.numeric(( Y %*% R.sample %*% Y) / mu2 )
                    EQ.sample <- sum(diag(R.sample))
                    trRR <- sum(R.sample*R.sample)
                    tr2R <- EQ.sample*EQ.sample
                    trR2 <- sum(diag(R.sample*R.sample))
                    if (model == 'logistic'){
                        K <- ( 1 - 6 * mu + 6 * mu * mu ) / mu2
                        varQ <- K * ( trR2 - tr2R / n ) + 2 * trRR - 2 * tr2R / (n-1)
                    }else{
                        varQ <- (2 / (n+1)) * ( (n-1) * trRR - tr2R )
                    }    
                    scl <- varQ / (2 * EQ.sample)
                    dfr <- EQ.sample / scl
                    ps[sample] <- pf ( (scl * dfr / Q.sample), 10^6, dfr )
               }
               sampledist[[m]]<-ps
            }
            comparative.p <- mean(sampledist[[m]] < rep(p.value, times=ndraws))
        }else 
            comparative.p <- NA
        

        # the returns for this geneset
        res[index, "p.val"] <- p.value
        res[index, "Q"] <- Q
        res[index, "EQ"] <- EQ
        res[index, "seQ"] <- seQ
        res[index, "comp.p"] <- comparative.p
      }else{
        mu2 <- NA
        res[index, "p.val"] <- NA
        res[index, "Q"] <- NA
        res[index, "EQ"] <- NA
        res[index, "seQ"] <- NA
        res[index, "comp.p"] <- NA
      }
  }
  
  # calculate influence per gene and expected influence
  XY <- t(X) %*% Y  
  obs.inf <- as.vector( sign(XY) * XY^2 / mu2 )
  exp.inf <- colSums(X * X)
  influence <- cbind(obs.inf, exp.inf)
  colnames(influence) <- c("observed","expected")

  # returns an object of type "gt.result"
  names(Y) <- samplenames
  new("gt.result",
    res = res, 
    X = X,
    Y = Y,
    influence = influence,
    test.genes = test.genes.lst)
}
#==========================================================


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
    if (gt@res[geneset, "test.n"] == 0)
      stop("empty pathway", call. = FALSE)
            
    # Recreate Y and R
    X <- gt@X[,test.genes]
    if (is.na(X))
      stop("empty pathway", call. = FALSE)
    m <- dim(X)[2]
    R <- (X %*% t(X)) / m
    Y <- gt@Y
    n <- length(Y)
    mu2 <- var(Y)
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
    for (i in 1:nchunks){
      Qs[(chunk * (i-1) + 1):(chunk * i)] <- permQ( R/mu2, Y, n, chunk ) 
      NULL
    }
    if (rest > 0)
      Qs[(chunk * nchunks + 1):nperm] <- permQ( R/mu2, Y, n, rest ) 

    # Draw histogram
    hst <- hist(Qs, xlim = c(0, 1.1 * max( c( Qs, Q ) ) ), 
      main = paste( "Histogram of Q for", nperm, "permutations of Y" ),
      xlab = "Values of Q for permuted Y" )
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

geneplot <- function(gt, geneset = NULL, genesubset = NULL, ...)
{
    # check correct input of gt
    if ( !is(gt, "gt.result"))
        stop("permutations should be applied to a globaltest result", call. = FALSE)
        
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
        if ( length(intersect(test.genes, genesubset)) > 0 ) {
            test.genes <- intersect(test.genes, gene.subset)
        }else{
            stop("genesubset is not a subset of the selected geneset")
        }
    }
    if (gt@res[geneset, "test.n"] == 0)
      stop("empty pathway", call. = FALSE)
    
    # Extract the calculated influence
    influence <- gt@influence[test.genes, 1]
    up <- (sign(influence) == 1)
    influence <- abs(influence)
    Einf <- gt@influence[test.genes, 2]
    m <- length(influence)
      
    # Plot reference line and influence per gene and color for up/down regulation
    plot( 0, xlim = c(1/2, m+1/2), ylim = c(0, 1.2 * max(influence, Einf) ), col = 0, xlab = "genenr", ylab = "influence" ,...)
    if (m <= 250) {
        rect(xleft = 1:m - 0.4, xright = 1:m + 0.4, ybottom = rep(0,times=m), ytop = influence, col = (up+2), border=0 )
    }else{
        lines(1:m, influence, lwd = 600 / m, type = 'h', col = (up+2))
    }
    lines( (1:(m+1))-1/2, c(Einf, Einf[m]), type = "s" )
    legend(1, 1.2 * max(c(influence,Einf)),
        c("pos. coregulated with Y", "neg. coregulated with Y"), fil = c(3,2))
        
    # Output: reference vector with gene labels if available
    invisible(names(influence))
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
    if ( !is(gt, "gt.result"))
        stop("permutations should be applied to a globaltest result", call. = FALSE)
        
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

    # recreate matrix R and find S = Y %o% Y
    X <- gt@X[,test.genes]
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
    plot(Sall, Rall, xlab = "Covariance between outcomes", ylab = "Covariance between expression patterns", col = 0,...)
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

checkerboard <- function(gt, geneset = NULL, sort = TRUE,...)
{   
    # check correct input of gt
    if ( !is(gt, "gt.result"))
        stop("permutations should be applied to a globaltest result", call. = FALSE)
        
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
    if (!is.logical(sort))
        stop("The option 'sort' should be either TRUE or FALSE", call. = FALSE)

    # recreate matrix R
    X <- gt@X[,test.genes]
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
    image(x = 1:n, y = 1:n, z = R>med, col = rainbow(2, v = c(0,1), s = c(1,0) ), xlab = label, ylab = label, 
        lab = c(n,n,50/n),...)
    par(pty = "s")
    invisible(newsamplenrs)
}
#==========================================================

#==========================================================
# Two functions to extract relevant information from 
# a gt.result object
#==========================================================
result <- function(gt)
{
    # check correct input of gt
    if ( !is(gt, "gt.result"))
        stop("results should be applied to a globaltest result", call. = FALSE)
    
    res <- gt@res
 
    if ( all(is.na(res[,"comp.p"])) ) {
      res <- res[,c("path.n","test.n","Q","EQ","seQ","p.val")]
    }else{
      res <- res[,c("path.n","test.n","Q","EQ","seQ","p.val","comp.p")]
    }
    
    res
}

p.value <- function(gt)
{
    # check correct input of gt
    if ( !is(gt, "gt.result"))
        stop("p.values should be applied to a globaltest result", call. = FALSE)

    gt@res[,"p.val"]
}

#==========================================================
setClass("gt.result", representation(
    res = "matrix", 
    X = "matrix",
    Y = "vector",
    influence = "matrix",
    test.genes = "list")
)


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

    #==========================================================
    # Function "show" prints out a result of type "gt.result"
    # such as results from a call to "globaltest"
    #==========================================================

    setMethod("show", "gt.result", function(object)
    {
      gt <- object
      npathways <- length(gt@test.genes)
      nsamples <- length(gt@Y)
      ngenes <- dim(gt@X)[2]

      cat("Global Test result:\n")
      cat("Data:", nsamples, "samples with", ngenes, "genes;", npathways, "pathways tested\n\n")
        
      res <- data.frame(gt@res)
      if ( all(is.na(res[,"comp.p"])) ) {
        res <- res[,1:6]
        colnames(res) <- c("genes","tested","Statistic Q","Expected Q","sd of Q","p-value")
      }else{
        colnames(res) <- c("genes","tested","Statistic Q","Expected Q","sd of Q","p-value","comp. p")
      }
      print(signif(res, digits = 5))
    }, where = where)
    #==========================================================


    invisible(NULL)
}
#==========================================================
