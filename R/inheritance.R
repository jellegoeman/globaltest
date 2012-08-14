inheritance <- function(test, sets, weights, ancestors, offspring, Shaffer, homogeneous=TRUE, trace) {
                                              
  # input checking 1: sets
  if (missing(sets) && is(test, "gt.object") && !is.null(test@subsets)) 
    sets <- test@subsets
  if (missing(sets)) stop("argument \"sets\" is missing with no default")
  if (is.character(sets))
    if (is(test, "gt.object") && !is.null(test@subsets)) {
      test <- test[sets]
      sets <- test@subsets
    }
  if(is(sets,"hclust")) 
    sets=as.dendrogram(sets)
  if(is(sets,"dendrogram")) {
    sets <- dendro2sets(sets)
    ancestors <- sets$ancestors
    sets <- sets$sets
  }
  if (is.null(names(sets)))
    stop("sets input has no names attribute.")
                                        
  # input checking 2: ancestors and offspring
  if (missing(ancestors) && is(test, "gt.object") && !is.null(test@structure$ancestors))
    ancestors <- test@structure$ancestors
  if (missing(offspring) && is(test, "gt.object") && !is.null(test@structure$offspring))
    offspring <- test@structure$offspring
  if (missing(ancestors) && missing(offspring)) {     # Infer from sets
    ancestors <- new.env(hash=TRUE)
    offspring <- new.env(hash=TRUE)
    for (i in 1:length(sets)) {
      truenamei <- names(sets)[i]
      namei <- as.character(i)      # workaround: variable names in environments may not be over 256 characters
      for (j in 1:length(sets)) {
        truenamej <- names(sets)[j]
        namej <- as.character(j)
        if (i != j && length(sets[[i]]) <= length(sets[[j]]) && all(sets[[i]] %in% sets[[j]])) {
          ancestors[[namei]] <- c(ancestors[[namei]], truenamej)
          offspring[[namej]] <- c(offspring[[namej]], truenamei)
        }
      }
    }
    ancestors <- as.list(ancestors)
    names(ancestors) <- names(sets)[as.numeric(names(ancestors))]
    offspring <- as.list(offspring)
    names(offspring) <- names(sets)[as.numeric(names(offspring))]
  }      
  if ((!missing(ancestors)) && is.environment(ancestors))
    ancestors <- as.list(ancestors)
  if ((!missing(offspring)) && is.environment(offspring))
    offspring <- as.list(offspring)
  if (missing(ancestors))
    ancestors <- turnListAround(offspring)
  if (missing(offspring))
    offspring <- turnListAround(ancestors) 
                               
  # get parents and children from ancestors and offspring
  parents <- ancestors2parents(ancestors)
  children <- ancestors2parents(offspring)
  leaves <- names(sets)[sapply(offspring[names(sets)], length) == 0]
                                        
  # prevent confusion between "weights" method and argument
  if (missing(weights)) weights <- NULL               
  if (is(weights, "gt.object") && is.list(test@weights) && length(test@weights) > 1)
    stop("The inheritance method is not applicable with individually weighted sets")
                              
  # Check whether Shaffer may be used
  if (is(test,"gt.object") && missing(Shaffer)) {
    shaffer.node <- names(sets)[sapply(names(sets), function(st) (!(st %in% leaves)) && all(children[[st]] %in% leaves))]
    Shaffer <- (length( setdiff(unlist(sets[shaffer.node]), unlist(sets[leaves])) ) == 0)
  }  
                                            
  # Check default trace
  if (missing(trace))    
    trace <- gt.options()$trace
                                   
  # Perform test for sets
  K <- length(sets)
  digitsK <- trunc(log10(K))+1
  if (is(test, "gt.object")) {
    # find out which tests have already been performed
    # sort subsets for quicker checking
    # for some reason a double for-loop is much quicker than a double sapply here!
    which.sort <- which(sapply(sets, length) %in% size(test))
    sets[which.sort] <- lapply(sets[which.sort], sort)
    test@subsets <- lapply(test@subsets, sort)
    found <- rep(NA, length(test))
    for (i in seq_along(test@subsets)) {
      subs <- test@subsets[[i]]
      for (j in seq_along(which.sort)) {
        jj <- which.sort[j]
        set <- sets[[which.sort[j]]]
        if ((length(subs) == length(set)) && all(subs == set)) {
          found[i] <- jj
          which.sort <- which.sort[-j]
          break
        } 
      }
    }
    # get test weights if available
    # NOTE: if inheritance was called from covariates(), test@weights is a list of length 1 
    # regardless of length(gt.object). We retrieve these test weigths and repair the object
    # TODO: program this in a more elegant way
    test.weights <- test@weights[[1]]
    # update object
    newgt <- test
    where.found <- match(1:K, found)
    progress <- sum(!is.na(where.found))
    newgt@result <- matrix(0,K,5)
    for (i in 1:K) {
      if (is.na(where.found[i])) {
        if (trace) {
          if (progress > 1)
            cat(rep("\b", digitsK+trunc(log10(progress))+4), sep="")
          progress <- progress +1
          cat(progress, "/", K)
          flush.console()
        }
        if (is.null(test.weights))
          newgt@result[i,] <- test@functions$test(sets[[i]])
        else
          newgt@result[i,] <- test@functions$test(sets[[i]], test.weights[sets[[i]]])
      } else {
        newgt@result[i,] <- test@result[where.found[i],]  
      }
    }
    rownames(newgt@result)<- names(sets)
    colnames(newgt@result) <- c("p-value", "Statistic", "Expected", "Std.dev", "#Cov")
    newgt@subsets <- sets
    newgt@extra <- NULL
    if (is.null(alias(test)) && (!all(names(newgt)[!is.na(where.found)] == names(test)[where.found[!is.na(where.found)]])))
      alias(test) <- names(test)
    if (!is.null(alias(test))) {
      alias(newgt) <- rep("", length(newgt))
      alias(newgt)[!is.na(where.found)] <- alias(test)[where.found[!is.na(where.found)]]
    }
    if (!is.null(test.weights))
      newgt@weights <- sapply(sets, function(st) test.weights[st])
    rawp <- p.value(newgt)
  } else {
    newgt <- NULL
    rawp <- sapply(1:K, function(i) {
      if (trace) {
        if (i > 1) cat(rep("\b", digitsK+trunc(log10(i-1))+4), sep="")
        cat(i, "/", K)
        flush.console()
      }
      test(sets[[i]])
    })
    names(rawp) <- names(sets)
  }

  # Get weights for inheritance: 1) get leaf weights if missing
  if (is.null(weights) && is(test, "gt.object")) {
    # most likely case: single top: first check
    tops <- which(sapply(ancestors[names(sets)], length) == 0)
    if (length(tops) == 1)
      weights <- weights(newgt[tops])
    else {
      temp <- test
      temp@subsets <- list(unique(unlist(sets[tops])))
      weights <- weights(temp)
    }
  } 
                                    
  # Get weights for inheritance: 2) get set weights
  if (is.null(weights)) {
    set.weights <- sapply(sets, length)
  } else {
    set.weights <- sapply(sets, function(st) sum(weights[st]))
  }

  #now call core function:
  adjustedP <- .inherit(ps=rawp, parents=parents, children=children, offspring=offspring, weights=set.weights, Shaffer=Shaffer, homogeneous =homogeneous)             # named vector

  if (trace==1) cat("\n")
  if (!is.null(newgt)) {
    extra <- newgt@extra
    extra[["inheritance"]] <- adjustedP
    newgt@extra <- as.data.frame(extra) 
    rownames(newgt@extra) <- names(newgt)
    newgt@structure <- list(ancestors=ancestors,offspring=offspring) 
    return(newgt)
  } else {
    return(data.frame(raw.p = rawp, inheritance = adjustedP))
  }

}

#############################################################################################################
################################################# FORMAT CONVERTING FUNCTIONS, SOME OTHER USEFUL FUNCTIONS and the core function


################################################################################
dendro2sets <- function(hc){
  ####################
  #create set names from the dendrogram
  do.sets <- function(x, parent.label, parent.branch, struct) {
    for (i in 1:length(x)) {
      parentnum <- as.character(which(struct$names == parent.label))
      currentnum <- as.character(length(struct$names)+1)
      newset <- labels(x[[i]])
      newname <- paste(parent.label,"[",i,sep="")
      if (is.leaf(x[[i]]))
        newname <- paste(newname, labels(x[[i]]), sep=":")
      newancestors <- c(struct$ancestors[[parentnum]], parent.label)
      newbranch <- c(parent.branch, i)
      struct$sets[[currentnum]] <- newset
      struct$ancestors[[currentnum]] <- newancestors
      struct$branch[[currentnum]] <- newbranch
      struct$names <- c(struct$names, newname)
      if (!is.leaf(x[[i]])) {
        struct <- do.sets(x=x[[i]], parent.label=newname, parent.branch=newbranch, struct)
      }
    }
    return(struct)
  }
  ####################

  struct <- list()
  struct$sets=new.env(hash=TRUE)
  struct$sets$"1" <- labels(hc)
  struct$ancestors <- new.env(hash=TRUE)
  struct$branch <- new.env(hash=TRUE)
  struct$branch$"1" <- numeric(0)
  struct$names <- "O"
  
  if (!is.leaf(hc))
    struct <- do.sets(hc, "O", numeric(0), struct)
  else
    struct$names <- paste(struct$names, labels(hc), sep=":")
  
  #convert to lists
  struct$sets <- as.list(struct$sets)
  names(struct$sets) <- struct$names[as.numeric(names(struct$sets))]
  struct$ancestors <- as.list(struct$ancestors)
  names(struct$ancestors) <- struct$names[as.numeric(names(struct$ancestors))]
  struct$branch <- as.list(struct$branch)
  names(struct$branch) <- struct$names[as.numeric(names(struct$branch))]
  
  return(struct)
}


#################################################
# returns the parent of element id
#get.parent<-function(id){
#  temp=sapply(structure$offspring[structure$ancestor[[id]]],length)
#  structure$ancestor[[id]][temp == min(c(unlist(temp),Inf))]
#  }
#################################################
# outputs the leaves of the structure
leaf.list<-function(structure) setdiff(names(structure$ancestor),names(structure$offspring))

#test if each element of all.list is a leaf node in the structure
is.leaf.list<-function(structure,all.list) {
  ans=rep(FALSE,length(all.list))
  names(ans)=all.list
  ans[intersect(leaf.list(structure) ,names(ans))]=TRUE
  ans
  }
  

######################################################################################################################
## here start the main function
.inherit <- function( ps,             # named vector
                      weights,        # named vector
                      parents,        # named list
                      children,
                      offspring,
                      Shaffer = TRUE,
                      homogeneous = FALSE) {

  # how many tests? 
  m <- length(ps)
  nms <- names(ps)
  weights <- weights[nms]
                                        
  # find top and leaves
  top <- sapply(parents[nms], length) == 0
  names(top) <- nms
  leaf <- sapply(children[nms], length) == 0
  names(leaf) <- nms
                          
  # find parents of leaves
  leaf.parents <- children[sapply(children, function(ch) any(ch %in% nms[leaf]))]
  leaf.siblings <- lapply(leaf.parents, function(lp) {
    nleaves <- sum(leaf[lp])
    if (nleaves > 1) lp else lp[!leaf[lp]]
  })

  # initialize
  basealpha <- rep(0,m)
  names(basealpha) <- nms
  basealpha[top] <- weights[top]/sum(weights[top])   # basealpha should add up to 1. Is multiplied by alpha later
  shaffer <- rep(1, m)
  rejected <- rep(FALSE, m)
  names(rejected) <- nms
  extinct <- weights == 0    # nodes with weight zero are never inherited to
  names(extinct) <- nms
  adjp <- rep(1, m)
  names(adjp) <- nms
  
  # a flag to do plain Meinshausen
  Meinshausen <- FALSE
  
  # start the procedure
  alpha <- 0
  ready <- FALSE
  while (!ready) {
                      
    # phase 1: reject
    testalpha <- basealpha * shaffer
    newly.rejected <- (ps/testalpha <= alpha) & (!rejected)
    newly.rejected[ps == 0 & testalpha == 0] <- FALSE     # do not reject when 0/0

    if (any(newly.rejected)) {
      adjp[newly.rejected] <- alpha
      rejected <- rejected | newly.rejected

      # phase 2: recalculate extinctness
      extinct[(!extinct) & rejected] <- sapply((1:m)[(!extinct) & rejected], function(i) 
        leaf[i] || all(rejected[offspring[[nms[i]]]]) || sum(weights[offspring[[nms[[i]]]]]) == 0)

      # phase 3: recalculate Shaffer
      if (Shaffer) {
        shaffer <- rep(1, m)
        names(shaffer) <- names(ps)
        if (homogeneous) {
          bottom <- setdiff(nms[rejected], unlist(parents[nms[rejected]]))
          bottom <- setdiff(bottom, nms[leaf])
          minweights <- unlist(lapply(bottom, function(bt) {
            min(weights[intersect(offspring[[bt]], nms[leaf])])
          }))
          sumweights <- sum(weights[leaf & !rejected])
          shaffer[1:m] <- sumweights / (sumweights - sum(minweights))
          for (bt in bottom) {
            if (length(children[[bt]]) > 1) {
              mw <- unlist(lapply(children[[bt]], function(ch) {
                min(weights[intersect(c(ch,offspring[[ch]]), nms[leaf])])
              }))  
              names(mw) <- children[[bt]]
              smallest <- names(which.min(mw))
              second <- min(mw[names(mw) != smallest])
              shaffer[smallest] <- sumweights / (sumweights - sum(minweights) + mw[smallest] - second)
            } else
              shaffer[children[[bt]]] <- Inf
          }
        } else {
          for (lp in names(leaf.parents)[rejected[names(leaf.parents)]]) {
            desc <- leaf.parents[[lp]]
            lsibs <- leaf.siblings[[lp]]
            if (!any(rejected[desc])) {
              lvs <- desc[leaf[desc]]
              if (length(desc) == 1)
                shaffer[desc] <- Inf
              else {             
                smallest <- names(which.min(weights[lvs])) 
                shaffer[lsibs] <- sum(weights[desc]) /  (sum(weights[desc]) - weights[smallest])
                if (smallest %in% lsibs) {
                  second <- min(weights[setdiff(lvs, smallest)])
                  shaffer[smallest] <- sum(weights[desc]) / (sum(weights[desc]) - second)
                }  
              }
            }
          }
        }
      }
                        
      # phase 4: inherit
      while (any(rejected & (basealpha > 0))) {
        for (i in 1:m) {
          if (rejected[i] && basealpha[i] > 0) {
                     
            # find heirs
            if (extinct[i])
              if (top[i])         # whole tree is rejected
                heirs <- nms[top & (!extinct)]
              else                # heirs of leaf nodes or nodes under which all leaves are rejected
                if (Meinshausen)
                  heirs <- character(0)
                else if (homogeneous) {
                  heirs <- nms[(basealpha > 0) & (1:m != i)]
                }  
                else
                  heirs <- parents[[nms[i]]]
            else {                 # heirs of internal nodes
              heirs <- children[[nms[i]]]
              heirs <- heirs[!extinct[heirs]]
            }  
                                         
            # inherit basealpha to heirs
            basealpha[heirs] <- basealpha[heirs] + basealpha[i] * weights[heirs]/sum(weights[heirs])
            basealpha[i] <- 0
          }
        }
      }
    } else {        # no new rejections: increase alpha
      if (all(rejected))
        ready <- TRUE
      else {
        next.rejection <- ps / (basealpha * shaffer)
        next.rejection[ps == 0 & basealpha == 0] <- Inf
        alpha <- min(next.rejection)
        ready <- (alpha >= 1)
        newalpha <- TRUE
      }
    }
  }
  adjp[!rejected] <- 1
  return(adjp)
}