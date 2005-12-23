#==========================================================
# The class "gt.barplot" stores the result of a call to
# geneplot or sampleplot
#==========================================================
setClass("gt.barplot", 
  representation(
    res = "matrix", 
    labelsize = "numeric",
    drawlabels = "logical",
    legend = "vector",
    colourCode = "vector"
  )
)

#==========================================================
setMethod("result", "gt.barplot",
            function(object) 
{
    res <- object@res
    colouring <- factor(res[,4])
    levels(colouring) <- object@colourCode[2:1]
    res <- data.frame(res[,1:3], z.score = z.score(object), colouring)
    res
})

#==========================================================
# The names method for "gt.barplot" 
#==========================================================
setMethod("names", "gt.barplot", 
            function(x) 
{
  rownames(x@res)
})    

#==========================================================
# The names<- method for "gt.barplot" 
#==========================================================
setMethod("names<-", "gt.barplot", 
            function(x, value) 
{
  rownames(x@res) <- value
  x
})            


#==========================================================
# A z.score method for "gt.barplot" 
#==========================================================
setMethod("z.score", "gt.barplot", 
            function(x) 
{
    res <- x@res
    (res[,1] - res[,2]) / res[,3]
})      

#==========================================================
setMethod("[", "gt.barplot", 
            function(x, i, j,...,drop) 
{
  if (all(i %in% rownames(x@res)) | all(i %in% 1:nrow(x@res))
    | (is.logical(i) & (length(i) == nrow(x@res)))) {
      x@res <- x@res[i, ,drop = FALSE] 
      x
  } else {
    stop("invalid subscript")
  }
})  

#==========================================================
setMethod("length", "gt.barplot", 
            function(x) 
{
  nrow(result(x))
})            

#==========================================================
setMethod("scale", "gt.barplot", 
            function(x, center = FALSE, scale = TRUE) 
{
  if (!missing(center) && center) {
    x@res[,1] <- x@res[,1] - x@res[,2]
    x@res[,2] <- x@res[,2] - x@res[,2]
  }
  if (missing(scale) || scale) {
    norm <- x@res[,3]
    norm[norm == 0] <- 1
    x@res[,1] <- x@res[,1] / norm
    x@res[,2] <- x@res[,2] / norm
    x@res[,3] <- x@res[,3] / norm
  }
  x
})            

#==========================================================
setMethod("show", "gt.barplot", 
            function(object) 
{
  show(result(object))        
})

#==========================================================
setMethod("plot", "gt.barplot",
  function(x,y,...) {
    plot.gt.barplot <- function(x, genesubset, drawlabels, labelsize, addLegend = TRUE, ...) {
      if (!missing(genesubset))
        x <- x[genesubset]
      if (missing(drawlabels))
        drawlabels <- x@drawlabels
      if (missing(labelsize))
        labelsize <- x@labelsize
      influence <- x@res[,1]
      Einf <- x@res[,2]
      sd.inf <- x@res[,3]
      up <- x@res[,4]
      legend <- rownames(x@res)
      m <- length(influence)
      rangebars <- max(0,influence, Einf) - min(0,influence)
      minplot <- min(0,influence) 
      maxplot <- max(0,influence, Einf) + 0.2 * rangebars
      if (drawlabels & !is.null(legend)){
      # check for space in margin of plot
        plot.new()  
        labwidth <- max(strwidth(legend,"inches",labelsize))
        margins <- par("mai")
        par(new = TRUE, "mai" = c(max(margins[1], labwidth*1.3), margins[2:4]))
        plot( 0, xlim = c(1/2, m+1/2), ylim = c(minplot, maxplot), col = 0, xlab = "", ylab = "influence", xaxt="n",...)
        axis(1, 1:length(legend), legend, cex.axis = x@labelsize, las=2)
        par("mai"=margins)
      }else
        plot( 0, xlim = c(1/2, m+1/2), ylim = c(minplot, maxplot), col = 0, xlab = "nr", ylab = "influence", ...)
      # plot influence bars
      if (m <= 300) {
        rect(xleft = 1:m - 0.4, xright = 1:m + 0.4, ybottom = rep(0,times=m), ytop = influence, col = (up+2), border=0 )
      } else {
        lines(1:m, influence, lwd = 600 / m, type = 'h', col = (up+2))
      }
      # plot Einf reference line
      lines( (1:(m+1))-1/2, c(Einf, Einf[m]), type = "s" )
      # plot sd.inf marks
      nlines <- trunc((influence - Einf) / sd.inf)
      nlines[sd.inf < 10^-3] <- 0
      for (index in 1:max(abs(nlines))) {
        high <- Einf + sign(nlines) * apply(cbind(abs(index), abs(nlines)), 1, min) * sd.inf
        segments( 1:m - 0.4, high, 1:m + 0.4, high)
      }
      # write the legend
      if (addLegend) { 
        legend(1/2, maxplot, x@legend[1:2], fil = c(3,2))
      }  
      invisible(NULL)
    }
    plot.gt.barplot(x, ...)
  }
)
