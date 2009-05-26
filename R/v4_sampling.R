#==========================================================
# Function "sampling" compares the p.value(s) found with 
# p-values from randomly generated "pathways" of the same size
# as the tested pathway.
#==========================================================
sampling <- function(gt, geneset, ndraws = 10^3)
{
  # check correct input of gt
  if ( !is(gt, "gt.result"))
    stop("sampling should be applied to a globaltest result", call. = FALSE)
      
  # check correct input of ndraws
  if ( !(ndraws > 0) )
    stop("option ndraws should be a positive integer", call. = FALSE)

  # extract the right geneset vector
  if (!missing(geneset))
    gt <- gt[geneset]
   
  model <- .model(gt)
  adjusted <- .adjusted(gt)
  p <- .nGenes(gt)

  setsizes <- unique(sapply(gt@genesets, length))
  setsizes <- setsizes[!(setsizes %in% c(0,p))]
  names(setsizes) <- paste("size", setsizes, sep="")
  setsizes <- as.list(setsizes)
  setsizes <- setsizes[!names(setsizes) %in% names(gt@samplingZs)]

  maxchunksize <- 10^6

  sampling.list <- lapply(setsizes, function(m) {
    chunk <- trunc(maxchunksize / min(m,p-m))
    nchunks <- trunc(ndraws / chunk)
    rest <- ndraws - nchunks * chunk
    allchunks <- c(as.list(rep(chunk, nchunks)), rest)
    zs <- unlist(lapply(allchunks, function(aChunk) {
      if (m < p/2) {
        randomgenesets <- replicate(aChunk, sample(p, m), simplify = FALSE) 
      } else {
        randomgenesets <- replicate(aChunk, - sample(p, p-m), simplify = FALSE)
      }
      randomgt <- gt
      randomgt@genesets <- randomgenesets
      if (model == "linear") {
        res <- .linearglobaltest(randomgt)
      } else if (model == "logistic") {
        res <- .logisticglobaltest(randomgt)
      } else  if (model == "survival") {
        if (adjusted) {
          res <- .adjustedsurvivalglobaltest(randomgt)
        } else {
          res <- .unadjustedsurvivalglobaltest(randomgt)
        }
      } else if (model == "multinomial") {
        if (adjusted) {
          res <- .adjustedmultinomialglobaltest(randomgt)
        } else {
          res <- .unadjustedmultinomialglobaltest(randomgt)
        }
      }       
      (res[,1] - res[,2]) / res[,3]
    }))
    zs
  })
  gt@samplingZs <- c(gt@samplingZs, sampling.list)
  
  if (.method(gt) == 2)
    zscore <- z.score(gt)
  else {
    if (model == "linear") {
      res <- .linearglobaltest(gt)
    } else if (model == "logistic") {
      res <- .logisticglobaltest(gt)
    } else  if (model == "survival") {
      if (adjusted) {
        res <- .adjustedsurvivalglobaltest(gt)
      } else {
        res <- .unadjustedsurvivalglobaltest(gt)
      }
    } else if (model == "multinomial") {
      if (adjusted) {
        res <- .adjustedmultinomialglobaltest(gt)
      } else {
        res <- .unadjustedmultinomialglobaltest(gt)
      }
    }
    zscore <- (res[,1] - res[,2]) / res[,3]
  }

  nTested <- .nTested(gt)
  comp.p <- sapply(1:.nPathways(gt), function(ix) {
    m <- nTested[ix]
    if ((m == 0) || (m == p) )
      out <- NA
    else {
      zs <- gt@samplingZs[[paste("size", m, sep = "")]]
      out <- mean(zs > zscore[ix])
    }
    out
  })
  if ("Comparative p" %in% colnames(gt@res))
    gt@res[,"Comparative p"] <- comp.p
  else {
    gt@res <- cbind(gt@res, comp.p)
    colnames(gt@res)[ncol(gt@res)] <- "Comparative p"
  }
  
  gt
}

#==========================================================
