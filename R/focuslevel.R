focusLevel <- function(test, sets, focus, ancestors, offspring,  
  stop = 1, atoms = TRUE, trace) {

  # rename input
  useAtoms <- atoms
  if (stop <= 1) {
    maxalpha <- stop
    stopafter <- Inf
  } else {
    maxalpha <- 1
    stopafter <- stop
  }       
  if (missing(trace))
    trace <- gt.options()$trace
       
  # input checking 1: sets
  if (missing(sets) && is(test, "gt.object") && !is.null(test@subsets)) 
    sets <- test@subsets
  if (missing(sets)) stop("argument \"sets\" is missing with no default")
  if (is.character(sets))
    if (is(test, "gt.object") && !is.null(test@subsets)) {
      test <- test[sets]
      sets <- test@subsets
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
      namei <- names(sets)[i]
      for (j in 1:length(sets)) {
        namej <- names(sets)[j]
        if (i != j && length(sets[[i]]) <= length(sets[[j]]) && all(sets[[i]] %in% sets[[j]])) {
          ancestors[[namei]] <- c(ancestors[[namei]], namej)
          offspring[[namej]] <- c(offspring[[namej]], namei)
        }
      }
    }
  }      
  if ((!missing(ancestors)) && is.environment(ancestors))
    ancestors <- as.list(ancestors)
  if ((!missing(offspring)) && is.environment(offspring))
    offspring <- as.list(offspring)
  if (missing(ancestors))
    ancestors <- turnListAround(offspring)
  if (missing(offspring))
    offspring <- turnListAround(ancestors)    

  # Preprocess gt.object input and calculate focus level p-values
  if (trace) cat("Getting raw p-values: ")
  K <- length(sets)
  digitsK <- trunc(log10(K))+1
  if (is(test, "gt.object")) {
    if (is.list(test@weights))
      stop("The focus level method is not applicable with individually weighted sets")
    rawgt <- test
    test <- function(set) rawgt@functions$test(set)[1]
    if (all(names(sets) %in% names(rawgt))) {
      rawgt <- rawgt[names(sets)]
    } else {
      result <- t(sapply(1:K, function(i) {
        if (trace) {
          if (i > 1) cat(rep("\b", digitsK+trunc(log10(i-1))+4), sep="")
          cat(i, " / ", K, sep="")
          flush.console()
        }
        rawgt@functions$test(sets[[i]])
      }))
      rawgt@result <- result
      colnames(rawgt@result) <- c("p-value", "Statistic", "Expected", "Std.dev", "#Cov")
      rownames(rawgt@result) <- names(sets)
      rawgt@subsets <- sets
    }
    rawp <- p.value(rawgt)[focus]
  } else {
    rawgt <- NULL
    rawp <- sapply(1:K, function(i) {
      if (trace) {
        cat(rep("\b", 2*digitsK+3), i, " / ", K, sep="")
        flush.console()
      }
      test(sets[[i]])
    })
    names(rawp) <- names(sets)
  }
  if (trace) cat(rep("\b", 2*digitsK+25), sep="")

     
  # Prepare the objects that store the test results
  pFocus <- rawp[focus]                                   # Stores the p-values of focus nodes
  sigFocus <- logical(length(focus))                # Keeps track of focus nodes which are declared significant
  names(sigFocus) <- focus
  emptyFocus <- logical(length(focus))              # Keeps track of subgraphs which are completely significant
  names(emptyFocus) <- focus
  atoms <- new.env(hash=TRUE, size=length(focus))   # Stores the atoms of offspring sets of each focus node used to make unions
  unions <- new.env(hash=TRUE, size=length(focus))  # Keeps track of all unions that may be tested
  #pUnions <- vector("list", length(focus))          # Stores the p-values of unions that may be tested
  #names(pUnions) <- focus 
  pUnions <- new.env(hash=TRUE, size=length(focus))
  
  sigUnions <- vector("list", length(focus))        # Keeps track of unions that have so far been called significant
  names(sigUnions) <- focus 
  offspringAtoms <- lapply(as.list(focus), function(term) { # Stores the offspring nodes as unions of atoms
    offnames <- intersect(offspring[[term]], names(sets))
    off <- matrix(,length(offnames),0)
    rownames(off) <- offnames
    off
  })
  names(offspringAtoms) <- focus
  sigOffspring <- vector("list", length(focus))     # Stores which offspring has already been declared significant
  names(sigOffspring) <- focus
  adjustedP <- rep(NA, length(sets))             # Prepare the vector of corrected p-values
  names(adjustedP) <- names(sets)

  # Find all GO terms above the focus level
  forefathers <- unique(unlist(ancestors[focus]))
  forefathers <- setdiff(forefathers, focus)
  forefathers <- setdiff(forefathers, unique(unlist(offspring[focus])))

  # Initialize
  alpha <- 0
  holm <- length(focus)
  indirectAffected <- list()
  ready <- FALSE
  change <- TRUE
  if (missing(trace)) trace <- interactive()

  
  while (! ready) {
  
    # Find the focus terms where the action is
    affected1 <- focus[(!sigFocus) & (pFocus * holm <= alpha)]    # Find newly significant focus level terms
    opensubgraphs <- focus[sigFocus & (!emptyFocus)]              # Find the subgraphs in which something may happen
    minP <- sapply(opensubgraphs, function(ff) min(pUnions[[ff]][!sigUnions[[ff]]]))  # Find the smallest p-value not yet declared significant
    if (length(minP) > 0) {
      names(minP) <- opensubgraphs
      affected2 <- names(minP * holm <= alpha)                      # Find those subgraphs for which that minimum is significant
    } else {
      affected2 <- character(0)
    }
    affected <- unique(c(affected1, affected2, names(indirectAffected)))  # Include those significant through upward propagation
    
    newSigGO <- character(0)
    change <- FALSE
    
    for (term in affected) {
      
      # If a focus term was declared significant. Start a new subgraph
      if (!sigFocus[term]) {                                
        sigFocus[term] <- TRUE
        newSigGO <- c(newSigGO, term)                       # A GO term has been declared significant
        if (trace > 1) cat("Rejected:", term, "\n")
        alloffspring <- intersect(offspring[[term]], names(sets))
        if (length(alloffspring) > 0) {
          offspringSets <- sets[alloffspring]
          if (useAtoms)
            TermAtoms <- getAtoms(offspringSets)                                  # Get the set of atoms
          else
            TermAtoms <- getNoAtoms(offspringSets) 
          atoms[[term]] <- TermAtoms
          unions[[term]] <- matrix(rep(TRUE, length(TermAtoms)), 1, length(TermAtoms))  # Prepare the matrix of unions of atoms
          sigUnions[[term]] <- FALSE
          pUnions[[term]] <- sapply(list(unique(unlist(TermAtoms))), test)                  # Start by testing the union of all atoms
          offspringAtoms[[term]] <- matrix(sapply(TermAtoms, function(y)            # Write all offspring as unions of atoms
            sapply(offspringSets, function(x) all(y %in% x))), ncol = length(TermAtoms))
          rownames(offspringAtoms[[term]]) <- names(offspringSets)
          sigOffspring[[term]] <- logical(length(offspringSets))
          names(sigOffspring[[term]]) <- names(offspringSets)
        } else {                                                                # An empty version for when the focus term is an end node
          sigUnions[[term]] <- logical(0)
          atoms[[term]] <- list()
          unions[[term]] <- matrix(,0,0)
          offspringAtoms[[term]] <- matrix(,0,0)
          sigOffspring[[term]] <- logical(0)
        }
        change <- TRUE
      } 
      
      # Propagate significance from offspring that was declared significant in another subgraph
      propagate <- indirectAffected[[term]]                 
      propagate <- propagate[!(propagate %in% term)]
      propagate <- propagate[!sigOffspring[[term]][propagate]]
      if (length(propagate) > 0) {
        pPatterns <- matrix(,0,length(atoms[[term]]))   # Get all superset patterns of a propagated GO term 
        for (i in 1:length(propagate)) {
          basePattern <- offspringAtoms[[term]][propagate[[i]],] # The propagated GO term as unions of atoms
          ancestorPatterns <- fillPattern(basePattern)      # All superset patterns as unions of atoms
          ancestorPatterns <- ancestorPatterns[!intersectPatterns(ancestorPatterns, pPatterns),,drop=FALSE]    # Remove duplicates
          pPatterns <- rbind(pPatterns, ancestorPatterns)
        }
        matchedPatterns1 <- intersectPatterns(pPatterns, unions[[term]])  # These new patterns were already testable
        matchedPatterns2 <- intersectPatterns(unions[[term]], pPatterns)  # These already testable patterns are now called significant
        newPatterns <- pPatterns[!matchedPatterns1,, drop=FALSE]          # These are new patterns
        if (trace > 1) cat("\t\tSignificance of", nrow(newPatterns), "sets from", propagate, " propagated to", term, "\n")
        unions[[term]] <- rbind(unions[[term]], newPatterns)              # Make all supersets testable and give them raw p 0
        pUnions[[term]][which(matchedPatterns2)] <- 0                     
        pUnions[[term]] <- c(pUnions[[term]], rep(0, nrow(newPatterns)))
        sigUnions[[term]] <- c(sigUnions[[term]], rep(FALSE, nrow(newPatterns)))    # These patterns will be called significant later
        change <- TRUE
      }
      
      # Find the testable non-significant unions with small p-values
      newsigs <- (pUnions[[term]] * holm <= alpha) & !sigUnions[[term]]     
      newsigs <- which(newsigs)
      
      # Find whether there are GO terms among the newly significant unions
      sigPatterns <- unions[[term]][newsigs,, drop=FALSE]
      offspringPatterns <- offspringAtoms[[term]][!sigOffspring[[term]],,drop=FALSE]     # Patterns of not already significant offspring
      matched <- intersectPatterns(offspringPatterns, sigPatterns)                  # Any new offspring terms among the significant patterns?
      newSigOffspring <- rownames(offspringPatterns)[matched]
      if (trace > 1 && (length(setdiff(newSigOffspring, propagate))>0)) {
        cat("Signficant:", setdiff(newSigOffspring, propagate), "\n")
      }
      newSigGO <- c(newSigGO, newSigOffspring)                 
      sigOffspring[[term]][newSigOffspring] <- TRUE
      
      # Expand the current subgraph by finding new unions that may now be tested, and find the p-values
      newpatterns <- matrix(, 0, length(atoms[[term]]))
      for (ix in newsigs) {
        sigUnionsMatrix <- unions[[term]][sigUnions[[term]], ,drop = FALSE]         # The unions so far significant for this term
        pattern <- unions[[term]][ix,]      # the significant pattern
        if (sum(pattern) > 1) {
          for (iy in (1:length(pattern))[pattern]) { # For-loop over the TRUEs
            newpattern <- pattern
            newpattern[iy] <- FALSE # All direct subsets have one extra FALSE
            reallyNew <- !any(intersectPatterns(matrix(newpattern,nrow=1), unions[[term]])) # Is the pattern really new?
            if (reallyNew) {      # Are the other parents of this pattern also present?
              parentspresent <- all(sapply((1:length(newpattern))[!pattern], function(iz) { # Loop over the FALSEs of pattern
                newpatternparent <- newpattern
                newpatternparent[iz] <- TRUE                        # The parents of newpattern have one extra TRUE
                any(apply(sigUnionsMatrix, 1, function(pttn) all(pttn == newpatternparent))) # Is this superset present among the significant tests?
              }))
              if (parentspresent) {
                newpatterns <- rbind(newpatterns, newpattern)
                rownames(newpatterns) <- NULL
              }
            }
          }
        }
        sigUnions[[term]][ix] <- TRUE # Only now call the term itself significant (prevents duplicate patterns)
      }
      unions[[term]] <- rbind(unions[[term]], newpatterns)
      sigUnions[[term]] <- c(sigUnions[[term]], rep(FALSE, nrow(newpatterns)))
      if (length(newpatterns) > 0) {                        # Calculate and store their p-values
        if (trace > 1) cat("\t\ttests:", nrow(newpatterns), "in", term, "\n")
        newsets <- lapply(as.list(1:nrow(newpatterns)), function(i) { # Assemble gene sets from the atoms
          unique(unlist(atoms[[term]][newpatterns[i,]]))
        })
        newpvalues <- sapply(newsets, test)
        pUnions[[term]] <- c(pUnions[[term]], newpvalues) 
        change <- TRUE
      }
      
      # Is the subgraph emptied?
      empty <- (length(atoms[[term]]) == 0) || (all(sigOffspring[[term]]))
      if (empty) {                                          
        emptyFocus[term] <- TRUE
      }
    } 
    
    if (length(newSigGO)>0) {
      change <- TRUE                # Don't change alpha yet; there may be propagation to be done
      adjustedP[newSigGO] <- alpha  # Adjusted P for newly affected GO terms
    }
    
    # Find all forefather terms that are now significant
    affectedGOforefathers <- unique(unlist(ancestors[newSigGO])) # Find GO terms above the focus level significant through propagation
    affectedGOforefathers <- intersect(affectedGOforefathers, forefathers)
    newlyAffectedGOforefathers <- intersect(affectedGOforefathers, names(adjustedP[is.na(adjustedP)]))
    adjustedP[newlyAffectedGOforefathers] <- alpha
    if (trace > 1 && (length(newlyAffectedGOforefathers) > 0))
      cat("Rejected because of rejected offspring:", newlyAffectedGOforefathers, "\n")
      
    # For all significant GO terms, see if they are also present in other subtrees  
    affectedGOwithParents <- lapply(as.list(newSigGO), function(nsg) { 
      present1 <- sapply(offspringAtoms, function(off) nsg %in% rownames(off))
      present2 <- focus %in% nsg
      focus[present1 | present2]
    })      # Creates a list of new significant GO terms, listing the subtrees they appear in
    names(affectedGOwithParents) <- newSigGO
    indirectAffected <- turnListAround(affectedGOwithParents)  # Reverses the list to a list of subtrees, listing the significant GO offspring terms
    
    # Output progress information
    if (trace == 1) {
      cat(paste(rep("\b", 48), collapse=""))
      cat("Alpha = ", format(alpha, digits=3, scientific=T, width= 10), ": rejected terms: ", format(sum(!is.na(adjustedP)), width=5), ".", sep="")
      flush.console()
    }
    
    # If nothing happened, increase alpha
    if (!change) { 
      allPs <- sort(c(pFocus, unlist(as.list(pUnions))))
      allPs <- allPs[allPs * holm > alpha]
      alpha <- allPs[1] * holm
      names(alpha) <- NULL
      if (trace > 1) cat("\tAlpha =", format(alpha, digits=3, scientific=T, width= 10), "\n")
    } else {
      if (trace > 1 && (holm > sum(!emptyFocus))) cat("Holm's factor:", sum(!emptyFocus), "\n")
      holm <- sum(!emptyFocus)
    }
    ready <- (holm == 0) || (alpha > maxalpha) || (sum(!is.na(adjustedP)) >= stopafter)
  } 
  
  # ready
  if (trace==1) {
    cat(paste(rep("\b", 48), collapse=""))
    cat("Alpha = ", format(1, digits=3, scientific=T, width= 10), ": rejected terms: ", format(length(adjustedP), width=5), ".\n", sep="")
  }
  if (alpha >= 1)
    adjustedP[is.na(adjustedP)] <- 1
  
  # return
  if (!is.null(rawgt)) {
    extra <- rawgt@extra
    extra[["focuslevel"]] <- adjustedP
    rawgt@extra <- as.data.frame(extra) 
    rawgt@structure <- list(ancestors=ancestors,offspring=offspring)
    return(rawgt)  
  } else {
    return(data.frame(raw.p = rawp, focuslevel = adjustedP))
  }
}


# Finds all ancestor patterns of a pattern
fillPattern <- function(pattern) {
  k <- sum(!pattern)
  out <- matrix(T, 2^k, length(pattern))
  if (k>0) {
    out[,!pattern] <- sapply(1:k, function(i) {
      rep(c(T,F), each = 2^(i-1), times = 2^(k-i))
    })
  }
  out
}


# Finds patterns in "patterns" which also occur in "matchto"
# Equivalent to %in% but for matrices
intersectPatterns <- function(patterns, matchto) {
  if (nrow(matchto) > 0) {
    matched <- apply(patterns, 1, function(pattern1) {
      any(apply(matchto, 1, function(pattern2) 
        all(pattern1 == pattern2)
      ))
    })
  } else {
    matched <- logical(nrow(patterns))
  }
  dim(matched) <- NULL
  matched
}


# Starts from a named list and "turns it around"
# Returns a named list of length(all elements of the list)
# Each element of the list gives the names of the original list elements
# that contained that element
turnListAround <- function(aList) {
  newlist <- new.env(hash=TRUE)
  objs <- names(aList)
  if (is.null(objs)) objs <- seq_along(alist)
  for (i in objs) {
    for (j in aList[[i]]) {
      newlist[[j]] <- c(newlist[[j]], i)
    }
  }
  as.list(newlist)
}
      


# This function extracts the smallest possible number of atoms from 
# sets, such that all sets can be built as unions of these atoms
# Note: These atoms lead to slightly conservative closed testing
# example: 1,2,23 is reduced to 1,2,3 and leads to unnessesary testing of 13 prior to testing 1.
# But: 1,12,13,23 is now reduced to 1,2,3 while otherwise it would not be reduced.
getAtoms <- function(allsets) {
          
  if (!is(allsets, "list")) {                                           # Coerce to a list
    allsets <- list(allsets)
  }
  allsets <- allsets[sapply(allsets, length) > 0]                       # Remove empty sets
  allsets <- lapply(allsets, function(set) set[!is.na(set)])            # Remove NAs
  setdiffs <- vector("list", length(allsets))                           # Prepare output
  i <- 0
  lastAtom <- character(0)
  while (length(allsets) > 0) {
    allsets <- lapply(allsets, function(set) setdiff(set, lastAtom))    # Remove probes in the last atom from every set
    allsets <- allsets[sapply(allsets, length) > 0]                     # Remove empty sets
    if (length(allsets) > 0) {
      i <- i+1
      allsets <- allsets[sort.list(sapply(allsets, length))]            # Sort by size
      lastAtom <- allsets[[1]]                                          # Next atom: smallest non-empty set
      setdiffs[i] <- allsets[1]                                         # Store the atom
    }
  } 

  setdiffs[1:i]
}

# No atoms variant: only throws away redundant sets
getNoAtoms <- function(sets) {
                     
  if (!is(sets, "list")) {                                     # Coerce to a list
    sets <- list(sets)
  }
  sets <- sets[sapply(sets, length) > 0]                       # Remove empty sets
  sets <- lapply(sets, function(set) set[!is.na(set)])         # Remove NAs
  K <- length(sets)
  redundant <- rep(FALSE, K)                                   # redundant sets are sets that are unions of other sets
  for (i in 1:K) {
    contained <- sapply(sets, function(set) all(set %in% sets[[i]]))
    redundant[i] <- all(sets[[i]] %in% unlist(sets[contained & (1:K != i) & !redundant]))
  }
  sets[!redundant]
}


findFocus <- function(sets, ancestors, offspring, maxsize = 10, atoms = TRUE) {

  useAtoms <- atoms

  if (is(sets, "gt.object") && !is.null(sets@subsets)) 
    sets <- sets@subsets
  if (missing(sets)) stop("argument \"sets\" is missing with no default")
  if (is.character(sets))
    if (is(test, "gt.object") && !is.null(test@subsets)) {
      test <- test[sets]
      sets <- test@subsets
    }
  if (missing(ancestors) && missing(offspring)) {     # Infer from sets
    ancestors <- new.env(hash=TRUE)
    offspring <- new.env(hash=TRUE)
    for (i in 1:length(sets)) {
      namei <- names(sets)[i]
      for (j in 1:length(sets)) {
        namej <- names(sets)[j]
        if (i != j && length(sets[[i]]) <= length(sets[[j]]) && all(sets[[i]] %in% sets[[j]])) {
          ancestors[[namei]] <- c(ancestors[[namei]], namej)
          offspring[[namej]] <- c(offspring[[namej]], namei)
        }
      }
    }
  }      
  if ((!missing(ancestors)) && is.environment(ancestors))
    ancestors <- as.list(ancestors)
  if ((!missing(offspring)) && is.environment(offspring))
    offspring <- as.list(offspring)
  if (missing(ancestors))
    ancestors <- turnListAround(offspring)
  if (missing(offspring))
    offspring <- turnListAround(ancestors)    
                        
  # Prune offspring
  offspring <- lapply(offspring, function(set) intersect(set, names(sets)))                      
                        
  # Find all nodes with at most 5*maxatoms offspring sets
  # and at most maxatoms atoms
  allnodes <- names(sets)
  doable <- sapply(offspring[allnodes], length) <= 5 * maxsize
  names(doable) <- names(sets)
  for (term in allnodes[doable]) {
    if (length(offspring[[term]]) > 0) {
      if (useAtoms)
        atoms <- getAtoms(sets[offspring[[term]]])
      else
        atoms <- getNoAtoms(sets[offspring[[term]]])
    } else {
      atoms <- character(0)
    }
    doable[term] <- (length(atoms) <= maxsize)
  }

  # Find all doable nodes with no doable ancestors
  alldoable <- allnodes[doable]
  ancestors <- ancestors[doable]
  redundant <- sapply(alldoable, function(id) {
    any(ancestors[[id]] %in% alldoable)
  })
  alldoable[!redundant]
}


# Finds the leaves of a significant graph
leaveNodes <- function(object, alpha=0.05, ancestors,type) {

  if (missing(ancestors)) {     # Infer from sets
    sets <- object@subsets
    ancestors <- new.env(hash=TRUE)
    for (i in 1:length(sets)) {
      namei <- names(sets)[i]
      for (j in 1:length(sets)) {
        namej <- names(sets)[j]
        if (i != j && length(sets[[i]]) <= length(sets[[j]]) && all(sets[[i]] %in% sets[[j]])) {
          ancestors[[namei]] <- c(ancestors[[namei]], namej)
        }
      }
    }
    ancestors <- as.list(ancestors)
  }

 if (missing(type)) type=names(object@extra)[1]
 
  focusP <- object@extra[[type]]      
  if (is.null(focusP)) stop("No focus level p-values")
  
  significant <- names(object)[focusP <= alpha]
  sign.anc <- unique(unlist(ancestors[significant]))
  significant <- setdiff(significant, sign.anc)
  
  object[significant]
}

# Takes sets to ancestors mapping; converts it to sets to parents mapping
ancestors2parents <- function(ancestors) {
  lapply(ancestors, function(anc) {
    setdiff(anc, unique(unlist(ancestors[anc])))
  })
}