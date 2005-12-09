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
  if (ncol(gt@PermQs) > 0)
    stop("sampling cannot be applied to a permutation globaltest result", call. = FALSE)
      
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
  setsizes <- setsizes[!names(setsizes) %in% names(gt@SamplingPs)]

  maxchunksize <- 10^6


  sampling.list <- lapply(setsizes, function(m) {
    chunk <- trunc(maxchunksize / min(m,p-m))
    nchunks <- trunc(ndraws / chunk)
    rest <- ndraws - nchunks * chunk
    allchunks <- c(as.list(rep(chunk, nchunks)), rest)
    ps <- unlist(lapply(allchunks, function(aChunk) {
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
        if (adjusted) {
          res <- .adjustedlogisticglobaltest(randomgt)
        } else {
          res <- .unadjustedlogisticglobaltest(randomgt)
        }
      } else  if (model == "survival") {
        if (adjusted) {
          res <- .adjustedsurvivalglobaltest(randomgt)
        } else {
          res <- .unadjustedsurvivalglobaltest(randomgt)
        }
      }
      res[,4]
    }))
    ps
  })
  gt@SamplingPs <- c(gt@SamplingPs, sampling.list)
  
  pvalue <- p.value(gt)
  nTested <- .nTested(gt)
  comp.p <- sapply(1:.nPathways(gt), function(ix) {
    m <- nTested[ix]
    if ((m == 0) || (m == p) )
      out <- NA
    else {
      ps <- gt@SamplingPs[[paste("size", m, sep = "")]]
      out <- mean(ps < pvalue[ix])
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
