#==========================================================
# Geneplot plots the influence of each gene on the outcome 
#   of the test statistic
# See help(geneplot) for details
#==========================================================

geneplot <- function(gt, geneset, genesubset, plot = TRUE, scale = FALSE, drawlabels = TRUE, labelsize = 0.6, ...)
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
  if ( !missing(genesubset)) {
    if ( is.character(genesubset) ) {
      test.genes <- intersect(match(genesubset, colnames(gt@eX)), geneset)
      if ( length(geneset) == 0 ) 
        stop("genesubset is not a subset of the selected geneset")
    }else{
      if (all(genesubset %in% 1:length(geneset)))
        geneset <- test.genes[genesubset]
      else
        stop("genesubset is not a subset of the selected geneset")
    }
  }

  # calculate influence per gene and expected influence
  genelist <- as.list(geneset)
  gt@genesets <- genelist
  model <- .model(gt)
  adjusted <- .adjusted(gt)
  if (model == "linear") {
    res <- .linearglobaltestgamma(gt)
  } else if (model == "logistic") {
    res <- .logisticglobaltest(gt)
  } else if (model == "multinomial") {
    if (adjusted) {
      res <- .adjustedmultinomialglobaltest(gt)
    } else {
      res <- .unadjustedmultinomialglobaltest(gt)
    }
  } else if (model == "survival") {
    if (adjusted) {
      res <- .adjustedsurvivalglobaltest(gt)
    } else {
      res <- .unadjustedsurvivalglobaltest(gt)
    }
  }
  if (model == "multinomial") {
    Y <- .Y(gt)
    Y <- Y / (rep(1, nrow(Y)) %o% sqrt(colSums(Y*Y)))
    covar <- gt@eX[geneset,] %*% Y 
    up <- apply(covar, 1, which.max)
  } else {
    up <- 2 - ((gt@eX[geneset,] %*% .Y(gt)) >= 0) 
  }
  res <- cbind(res[,1:3,drop = FALSE], up)
  colnames(res) <- c("Influence", "Expected", "SD", "UP")
  
 
  # Output: gt.barplot object
  if (model == 'linear') {
    nameY <- as.character(.formula(gt)[[2]])
    colourCode <- c("+", "-")
    colour <- c(3,2)
    if (adjusted) {
      legend <- c(paste("positive correlation with residual", nameY),
        paste("negative correlation with residual", nameY))
    } else {
      legend <- c(paste("positive correlation with", nameY),
        paste("negative correlation with", nameY))
    }
  } else if (model == 'logistic') {
    colourCode <- c(paste("high in", .levels(gt)[1]), paste("high in", .levels(gt)[2]) )
    colour <- c(3,2)
    legend <- c(paste("higher expression in", .levels(gt)[1], "samples"), 
      paste("higher expression in", .levels(gt)[2], "samples"))
  } else if (model == 'survival') {
    colourCode <- c("+", "-")
    legend <- c(paste("positively associated with survival"), 
      paste("negatively associated with survival"))
  } else if (model == "multinomial") {
    colourCode <- sapply(.levels(gt), function(level) paste("high in", level))
    colour <- 1 + 1:length(.levels(gt))
    legend <- sapply(.levels(gt), function(level) paste("highest expression in", level, "samples"))
  }

  gtbar <- new("gt.barplot",     
    res = res,
    labelsize = labelsize, 
    drawlabels = drawlabels,
    legend = legend,
    colour = colour,
    colourCode = colourCode)
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
