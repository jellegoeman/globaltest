#==========================================================
# The function regressionplot allows the evaluation of 
#   possibly outlying samples.
# See help(regressionplot) for details
#==========================================================

regressionplot <- function(gt, geneset, sampleid, ...)
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

    # recreate matrix R and find S = Y %o% Y
    X <- gt@eX[geneset,,drop=FALSE]
    Y <- .Y(gt)
    n <- .nSamples(gt)
    R <- as.vector(crossprod(X))
    if (.model(gt) == "multinomial") {
      S <- rep(0,n*n)
      for (ix in 1:ncol(Y))
        S <- S + as.vector(Y[,ix] %o% Y[,ix])
    } else
      S <- as.vector(Y %o% Y)
    
    
    # Check correct input of sampleid
    if ( missing(sampleid) ) {
      sampleid <- NULL
    }else{
      if (!( all(sampleid %in% 1:n) | all(sampleid %in% colnames(X)) ))
        stop("Option sampleid incorrect", call. = FALSE)
    }
            
    # Extract relevant entries from S and R
    samples <- rep(FALSE,times = n)
    names(samples) <- colnames(X)
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
      points(Srest, Rrest,col = 4, cex=0.5)
      abline(lm(Rall ~ Sall), col = 4)
    }
    if (length(Rsub) > 0){
      points(Ssub, Rsub, col = 2, pch = 4, cex = 2)
      abline(lm(Rsub ~ Ssub), col = 2, lty = 2)
    }
    
    # No output
    invisible(NULL)
}   
#==========================================================
