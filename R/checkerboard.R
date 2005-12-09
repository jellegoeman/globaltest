#==========================================================
# Function Checkerboard visualizes a result of globaltest
#   by visualizing the covariance matrix between the different samples
# See help(checkerboard) for details
#==========================================================

checkerboard <- function(gt, geneset, sort = TRUE, drawlabels = TRUE, labelsize = 0.6,...)
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
    if (!is.logical(sort))
        stop("The option 'sort' should be either TRUE or FALSE", call. = FALSE)

    # recreate matrix R
    X <- gt@eX[geneset,,drop=FALSE]
    Y <- .Y(gt)
    R <- crossprod(X)
        
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
    rownames(newsamplenrs) <- colnames(X)

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
