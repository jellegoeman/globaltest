require("methods")

#==========================================================
# Function "globaltest" performs the global test
# X is a data matrix (n rows are samples, p columns are genes).
# NB: if the dimension of X does not fit the dimension of Y,
#       but t(X) does, t(X) is used.
# Code missing values as NA
# Y is a vector of n clinical outcomes.

# OPTIONS:
# model = 'logistic': for two-valued Y (default).
# model = 'linear': for continuous Y.
# test.genes = vector: to only test a subset of the genes.
#       vector can be a length p vector of 0 and 1
#           (1 = selected; 0 = not selected).
#       or a vector with the numbers of the selected genes
#           e.g. c(5,8,9) to test genes 5, 8 and 9.
#       default: test all genes.
# permutation = FALSE: give asymptotic version of test (default)
# permutation = TRUE: give permutation version of test.
# nperm = number: number of permutations sampled (default = 10^4).
#==========================================================

globaltest <- function(X, Y,
                        test.genes = NULL,
                        model = 'logistic',
                        permutation = FALSE,
                        nperm = NULL)

{
    # check for correct input of X and Y:
    # 1: coerce X into a matrix
    if (is(X, "exprSet"))
        X <- exprs(X)
    if (!is.vector(X) & !is.data.frame(X) & !is.matrix(X))
        stop("X should be of type 'matrix'.", call. = FALSE)
    X <- as.matrix(X)

    # 2: coerce Y into a vector
    if ( !is.vector(Y) & !is.matrix(Y) & !is.data.frame(Y) )
        stop("Y should be of type 'vector'.", call. = FALSE)
    if ( is.data.frame(Y) )
        Y <- as.matrix(Y)
    if ( is.matrix(Y) ) {
        if ( dim(Y)[1] == 1 )
            Y <- t(Y)
        if ( dim(Y)[2] == 1 ) {
            samples <- rownames(Y)
            Y <- as.vector(Y)
            names(Y) <- samples
        } else
            stop("Y should be a single vector.", call. = FALSE)
    }

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
    samplenames <- rownames(X)
    if (  ( is.null(samplenames) ) & ( !is.null( names(Y) ) ) )
        samplenames <- names(Y)
    if ( ( !is.null(rownames(X) ) ) & ( !is.null( names(Y) ) ) & ( !all(names(Y) == rownames(X)) ) )
        warning("Sample names in X inconsistent with sample names in Y.", call. = FALSE)
    dimnames(X) <- NULL
    names(Y) <- NULL

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

    # check for correct input of test.genes
    # 1: set default
    if (is.null(test.genes))
        test.genes <- rep(TRUE,times = p)
    else{
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
            names(test.genes) <- NULL
        }
        else
            if ( all(test.genes == as.integer(test.genes)) & (min(test.genes) >= 1) & (max(test.genes) <= p) ) {
                falses <- rep(FALSE, times = p)
                falses[test.genes] <- TRUE
                if ( length(test.genes) == sum(falses) )
                    test.names <- names(sort(test.genes))
                else
                    test.names <- NULL
                test.genes <- falses
            }
            else
                stop("Option test.genes should have 0 or 1 for each gene\n or be a list of gene numbers.", call. = FALSE)

        # 4: check compatibility of names
        genenames <- genenames[test.genes]
        if ( ( !is.null(genenames ) ) & ( !is.null( test.names ) ) & ( !all( genenames == test.names) ) )
           warning("Gene names in X inconsistent with gene names in test.genes.", call. = FALSE)

    }

    # check for correct input of permutation and nperm
    if ( (!is.logical(permutation)) | (length(permutation) != 1) )
        stop("Option permutation should be either TRUE or FALSE", call. = FALSE)
    if (!is.null(nperm)){
        if ( !(nperm == as.integer(nperm)) | (nperm < 1) | (length(nperm) != 1) )
            stop("Option nperm should be a single positive integer", call. = FALSE)
        else
            permutation <- TRUE
    }else{
        if (permutation)
            nperm <- 10^4
    }

    # center Y
    mu <- mean(Y)
    Y <- Y - mu

    # select genes to be tested
    X <- as.matrix(X[,test.genes])
    m <- sum(test.genes)

    # Center genes and set missing values to zero
    col.meanX <- colMeans(X, na.rm = TRUE)
    X.cent <- X - rep(1,times=n) %o% col.meanX
    X.cent[is.na(X.cent)] <- 0

    # calculate test statistic Q
    R <- ( X.cent %*% t(X.cent) ) / m
    mu2 <- switch(model,
        'linear' = var(Y),
        'logistic' = mu * (1-mu))
    Q <- as.numeric( ( Y %*% R %*% Y ) / mu2 )

    # Expectation, variance and p-value of Q from asymptotic formula's
    if (!permutation){
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
        Qs <- numeric(0)
    }

    # Interior function to calculate Q for a number npm of permutations of Y
    permQ <- function(stR, YY, npm){
        pms <- apply( matrix(rnorm(n * npm), n, npm), 2, sort.list )
        YY.pm <-  matrix( YY[pms], n, npm )
        colSums(( stR %*% YY.pm ) * YY.pm)
    }

    # Expectation, variance and p-value of Q from permutation version
    if (permutation)
    {
        Qs <- numeric(nperm)
        chunk <- 5000
        nchunks <- trunc(nperm / chunk)
        rest <- nperm - nchunks * chunk
        for (i in 1:nchunks){
            Qs[(chunk * (i-1) + 1):(chunk * i)] <- permQ( R/mu2, Y, chunk )
            NULL
        }
        if (rest > 0)
            Qs[(chunk * nchunks + 1):nperm] <- permQ( X.cent, Y, rest ) / (m * mu2)
        EQ <- mean(Qs)
        seQ <- sd(Qs)
        p.value <- sum(Qs >= rep(Q, times=nperm)) / nperm
    }

    # calculate influence per gene and expected influence
    XY <- t(X.cent) %*% Y
    influence <- as.vector( sign(XY) * XY^2 / mu2 )
    exp.influence <- colSums(X.cent * X.cent)

    # returns an object of type "gt.result"
    names(influence) <- genenames
    names(Y) <- samplenames
    new("gt.result",
        Q = Q,
        EQ = EQ,
        seQ = seQ,
        p.value = p.value,
        Qs = Qs,
        matrixR = R,
        Y = Y,
        influence = influence,
        exp.influence = exp.influence,
        test.genes = test.genes)
}
#==========================================================


#==========================================================
# Function "permutations" gives a histogram of the values of Q
#   calculated for permutations of the clinical outcome
# An arrow points out the value of the true Q, for comparison
#==========================================================

permutations <- function(testresult)
{
    nperm <- length(testresult@Qs)
    if ( nperm > 0 ){
        hst <- hist(testresult@Qs, xlim = c(0, 1.1 * max( c( testresult@Qs, testresult@Q ) ) ),
            main = paste( "Histogram of Q for", nperm, "permutations of Y" ),
            xlab = "Values of Q for permuted Y" )
        h <- max(hst$counts)
        arrows( testresult@Q, h/5 , testresult@Q, 0 )
        text( testresult@Q, h/4.5, 'Q' )
        invisible(NULL)
    }else{
        stop("The function 'permutations' requires that\n the permutation version of 'globaltest' was used.", call. = FALSE)
    }
}
#==========================================================


#==========================================================
# Geneplot plots the influence of each gene on the outcome
#   of the test statistic
# See help(geneplot) for details
#==========================================================

geneplot <- function(testresult)
{
    if (!is(testresult, "gt.result"))
        stop("The function geneplot can only be applied to results of the function globaltest.", call. = FALSE)
    m <- length(testresult@influence)
    influence <- abs(testresult@influence)
    up <- ( sign(testresult@influence) == 1 )
    Einf <- testresult@exp.influence

    # Plot reference line and influence per gene and color for up/down regulation
    plot( 0, xlim = c(1/2, m+1/2), ylim = c(0, 1.2 * max(influence, Einf) ), col = 0, xlab = "genenr", ylab = "influence" )
    if (m <= 250)
        rect(xleft = 1:m - 0.4, xright = 1:m + 0.4, ybottom = rep(0,times=m), ytop = influence, col = (up+2), border=0 )
    else
        lines(1:m, influence, lwd = 600 / m, type = 'h', col = (up+2))
    lines( (1:(m+1))-1/2, c(Einf, Einf[m]), type = "s" )
    legend(1, 1.2 * max(c(influence,Einf)),
        c("pos. coregulated with Y", "neg. coregulated with Y"), fil = c(3,2))

    # Output: reference list with gene labels if available
    genenames <- names(influence)
    if (is.null(genenames))
        invisible(NULL)
    else
        if (m <= 50)
            data.frame(genenames)
        else
            invisible( data.frame(genenames) )
}
#==========================================================


#==========================================================
# The function regressionplot allows the evaluation of
#   possibly outlying samples.
# See help(regressionplot) for details
#==========================================================

regressionplot <- function(testresult, samplenr=NULL)
{
    if (!is(testresult, "gt.result"))
        stop("The function regressionplot can only be applied to results of the function globaltest.", call. = FALSE)

    # The matrices S and R to be plotted against each other
    R <- as.vector(testresult@matrixR)
    S <- as.vector(testresult@Y %o% testresult@Y)
    n <- length(testresult@Y)

    # Check correct input of samplenr
    if ( is.null(samplenr) )
        samplenr <- 1:n
    if ( !is.vector(samplenr) | !all(samplenr == as.integer(samplenr)) | (max(samplenr) > n) | (min(samplenr) < 1) )
        stop( paste("Option samplenr should be a vector of sample nrs beween 1 and ", n, sep = ""), call. = FALSE )

    # Extract relevant entries from S and R
    samples <- rep(FALSE,times = n)
    samples[samplenr] <- TRUE
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
    plot(Sall, Rall, xlab = "Covariance between outcomes", ylab = "Covariance between expression patterns", col = 0)
    if (length(Rrest) > 0){
        points(Srest, Rrest,col = 4)
        abline(lm(Rall ~ Sall), col = 4)
    }
    points(Ssub, Rsub, col = 2, pch = 4, cex = 1.5)
    abline(lm(Rsub ~ Ssub), col = 2)

    # Some explanation
    samplenames <- names( testresult@Y )
    if ( ( !all(samples) ) & ( !is.null(samplenames) ) )
        cat("Samples investigated:\n", samplenames[samplenr], "\n")

    # No output
    invisible(NULL)
}
#==========================================================


#==========================================================
# Function Checkerboard visualizes a result of globaltest
#   by visualizing the covariance matrix between the different samples
# See help(checkerboard) for details
#==========================================================

checkerboard <- function(testresult, sort = TRUE)
{
    if (!is(testresult, "gt.result"))
        stop("The function checkerboard can only be applied to results of the function globaltest.", call. = FALSE)
    if (!is.logical(sort))
        stop("The option 'sort' should be either TRUE or FALSE", call. = FALSE)

    # Sort Y if needed and return new samplenrs
    Y <- testresult@Y
    n <- length(Y)
    perm <- sort.list(Y)
    rperm <- sort.list(Y[n:1])
    if ( any(perm != 1:n ) & any(rperm != 1:n) & (sort)){
        cat("Samples in the checkerboard have been sorted by clinical outcome:\n")
        newsamplenrs <- t( matrix( c(1:n, sort.list( perm )), n, 2 ) )
        label <- "sorted samplenr"
        R <- testresult@matrixR[perm,perm]
    }else{
        sort = FALSE
        label <- "samplenr"
        R <- testresult@matrixR
        newsamplenrs <- t( matrix( c(1:n, 1:n), n, 2 ) )
    }
    rownames(newsamplenrs) <- c("samplnr.old", "samplenr.new")
    colnames(newsamplenrs) <- names(Y)

    # Calculate median non-diagonal element of R
    lowertriangle <- outer( 1:n, 1:n, ">" )
    med <- median(R[lowertriangle])

    # Draw plot
    image(x = 1:n, y = 1:n, z = R>med, col = rainbow(2, v = c(0,1), s = c(1,0) ), xlab = label, ylab = label,
        lab = c(n,n,50/n))
    par(pty = "s")
    if (!sort)
        invisible(newsamplenrs)
    else
        newsamplenrs
}
#==========================================================


#==========================================================
# gene2ix takes the genenames of the genes of interest
# and compares them to the genenames on the chip
# to get the indices of the genes of interest on the chip
#==========================================================

gene2ix <- function(genes, chip, logic = FALSE)
{
    g2x <- function(gns, chp)
    {
        lg <- length(gns)
        if ( lg < 100 )
            apply( outer(chp, gns, "=="), 1, any )
        else{
            half <- round(lg/2)
            ix1 <- g2x(gns[1:half], chp)
            ix2 <- g2x(gns[(half+1):lg], chp)
            ix1 | ix2
        }
    }

    if (!is.vector(chip) | !is.vector(genes))
        stop("Both chip and genes should be of type 'vector'.", Call.= FALSE)
    ix <- g2x( genes, chip )
    cat( sum(ix), "out of", length(genes), "genes present on the chip of", length(chip), "genes.\n" )
    if (logic) {
        names(ix) <- chip
        ix
    }
    else {
        allgenes <- 1:length(chip)
        names(allgenes) <- chip
        allgenes[ix]
    }
}
#==========================================================


#==========================================================
# .First.lib is called when the package is loaded
# It initializes the object "gt.result"
#   and the its "show" function.
#==========================================================

.First.lib <- function(libname, pkgname)
{
    setClass("gt.result", representation(
        Q = "numeric",
        EQ = "numeric",
        seQ = "numeric",
        p.value = "numeric",
        Qs = "vector",
        matrixR = "matrix",
        Y = "vector",
        influence = "vector",
        exp.influence = "vector",
        test.genes = "vector"))

    #==========================================================
    # Function "show" prints out a result of type "gt.result"
    # such as results from a call to "globaltest"
    #==========================================================

    setMethod("show", "gt.result", function(object)
    {
        # Print numbers of genes, samples and tested genes
        m <- length(object@influence)
        n <- length(object@Y)
        p <- length(object@test.genes)
        cat("Global Test result:\n")
        cat( m, "out of", p, "genes used;", n, "samples\n\n")

        # Print p-value
        p.value <- signif(object@p.value, digits = 4)
        cat( "p value =", p.value, "\n" )
        if (length(object@Qs) == 0)
            cat("    based on theoretical distribution\n\n")
        else
            cat("    based on", length(object@Qs), "permutations\n\n")

        # Print summary of test statistic
        Q <- signif(object@Q, digits = 4)
        EQ <- signif(object@EQ, digits = 4)
        seQ <- signif(object@seQ, digits = 4)
        cat( "Test statistic Q =", Q, "\n" )
        cat( "    with expectation EQ =", EQ, "\n" )
        cat( "    and standard deviation sdQ =", seQ, "under the null hypothesis\n" )
    })
    #==========================================================

    invisible(NULL)
}
#==========================================================
