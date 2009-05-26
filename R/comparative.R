#==========================================================
# Function "sampling" compares the p.value(s) found with
# p-values from randomly generated "pathways" of the same size
# as the tested pathway.
#==========================================================
comparative <- function(object, N = 1e3, z.scores = TRUE, trace)
{
  p <- object@functions$df()[3]

  if (is.null(object@weights)) {
  
    # comparative proportion can be calculated at once for all subsets with the same size
    sizes <- unique(size(object))
    sizes <- as.list(sizes)
  
    K <- length(sizes)
    digitsK <- trunc(log10(K))+1

    if (missing(trace)) trace <- interactive() && K > 1
  
    sampling.list <- lapply(seq_along(sizes), function(i) {
      if (trace) {
        cat(rep("\b", 2*digitsK+3), i, " / ", K, sep="")
        flush.console()
      }
      m <- sizes[[i]]
      if (m %in% c(0,p)) 
        pps <- NA
      else {      
        if (m < p/2) {
          randomsets <- replicate(N, sample(p, m), simplify = FALSE)
        } else {
          randomsets <- replicate(N, - sample(p, p-m), simplify = FALSE)
        }               
        pps <- sapply(1:N, function(j) {
          res <- object@functions$test(randomsets[[j]], calculateP = !z.scores)
          if (z.scores) 
            -(res["S"] - res["ES"]) / res["sdS"]     #minus z as proxy for p
          else
            res["p"]
        })
      }
    })
    names(sampling.list) <- paste("size", sizes, sep="")
  
    if (z.scores) true.p <- -z.score(object) else true.p <- p.value(object)
    true.size <- paste("size", size(object), sep="")
    comparative <- sapply(1:length(object), function(i) {
      if (true.size[i] %in% c(0,p))
        NA
      else
        mean(true.p[i] >= sampling.list[[true.size[i]]])
    })
  } else { 
    K <- length(object)
    digitsK <- trunc(log10(K))+1

    if (missing(trace)) trace <- interactive() && K > 1
  
    if (z.scores) true.p <- -z.score(object) else true.p <- p.value(object)
        
  
    if (is.null(object@subsets)) {
      comparative <- sapply(1:length(object), function(i) {
        randomweights <- replicate(N, object@weights[[i]][sample(p,p)], simplify=FALSE)
        if (trace) {
          cat(rep("\b", 2*digitsK+3), i, " / ", K, sep="")
          flush.console()
        }
        pps <- sapply(1:N, function(j) {
          if (trace) {
            cat(rep("\b", 2*digitsK+3), N*(i-1) + j, " / ", K, sep="")
            flush.console()
          }
          res <- object@functions$test(weight = randomweights[[j]], calculateP = !z.scores)            
          if (z.scores) 
            -(res["S"] - res["ES"]) / res["sdS"]     #minus z as proxy for p
          else
            res["p"]
        })
        mean(true.p[i] >= pps)
      })  
    } else {
      comparative <- sapply(1:length(object), function(i) {
        size <- length(object@subsets[[i]])
        randomsubsets <- replicate(N, sample(p,size), simplify=FALSE)
        if (trace) {
          cat(rep("\b", 2*digitsK+3), i, " / ", K, sep="")
          flush.console()
        }
        pps <- sapply(1:N, function(j) {
          if (trace) {
            cat(rep("\b", 2*digitsK+3), N*(i-1) + j, " / ", K, sep="")
            flush.console()
          }
          res <- object@functions$test(subset = randomsubsets[[j]], weight = object@weights[[i]], calculateP = !z.scores)
          if (z.scores) 
            -(res["S"] - res["ES"]) / res["sdS"]     #minus z as proxy for p
          else
            res["p"]
        })
        mean(true.p[i] >= pps)
      })  
    }      
  }
  if (trace) cat("\n")
 
  if (is.null(object@extra))
    object@extra <- as.data.frame(matrix(,length(object),0), row.names=names(object))
  object@extra[["comparative"]] <- comparative
  object

}


#==========================================================
