gtGO <- function(..., GO, focus, maxalpha = 0.05, stopafter = 100, verbose = FALSE) {

  # check input
  if (missing(GO) || !(is(GO, "GOstructure"))) stop("Provide GO as a GOstructure object")
  if (!all(focus %in% GO@ids)) stop("Not all focus terms found in the GOstructure object")

  # Get the raw p-values
  smallgt <- globaltest(...)
  if (smallgt@method > 3) {
    stop("gtGO is not available with permutation testing")
  }
  rawgt <- .globaltest(smallgt, genesets(GO)[focus])

  # Prepare the objects that store the test results
  pFocus <- rawgt[,4]                               # Stores the p-values of focus nodes
  sigFocus <- logical(length(focus))                # Keeps track of focus nodes which are declared significant
  names(sigFocus) <- focus
  emptyFocus <- logical(length(focus))              # Keeps track of subgraphs which are completely significant
  names(emptyFocus) <- focus
  atoms <- vector("list", length(focus))            # Stores the atoms of offspring sets of each focus node used to make unions
  names(atoms) <- focus
  unions <- vector("list", length(focus))           # Keeps track of all unions that may be tested
  names(unions) <- focus 
  pUnions <- vector("list", length(focus))          # Stores the p-values of unions that may be tested
  names(pUnions) <- focus 
  
  sigUnions <- vector("list", length(focus))        # Keeps track of unions that have so far been called significant
  names(sigUnions) <- focus 
  offspring <- lapply(as.list(focus), function(term) { # Stores the offspring nodes as unions of atoms
    offnames <- unlist(GO@offspring[term])
    off <- matrix(,length(offnames),0)
    rownames(off) <- offnames
    off
  })
  names(offspring) <- focus
  sigOffspring <- vector("list", length(focus))     # Stores which offspring has already been declared significant
  names(sigOffspring) <- focus
  adjustedP <- rep(NA, length(GO))             # Prepare the vector of corrected p-values
  names(adjustedP) <- GO@ids

  # Find all GO terms above the focus level
  forefathers <- unique(unlist(GO@ancestors[focus]))
  forefathers <- setdiff(forefathers, focus)
  forefathers <- setdiff(forefathers, unique(unlist(GO@offspring[focus])))

  # Initialize
  alpha <- 0
  holm <- length(focus)
  indirectAffected <- list()
  ready <- FALSE
  change <- TRUE
  
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
        if (verbose) cat("Signficant:", term, "\n")
        alloffspring <- GO@offspring[[term]]
        if (length(alloffspring) > 0) {
          offspringSets <- GO@genesets[alloffspring]
          TermAtoms <- getAtoms(offspringSets)                                  # Get the set of atoms
          atoms[[term]] <- TermAtoms
          unions[[term]] <- matrix(rep(TRUE, length(TermAtoms)), 1, length(TermAtoms))  # Prepare the matrix of unions of atoms
          sigUnions[[term]] <- FALSE
          pUnions[[term]] <- .globaltest(smallgt, list(unlist(TermAtoms)))[,4]
          offspring[[term]] <- matrix(sapply(TermAtoms, function(y)            # Write all offspring as unions of atoms
            sapply(offspringSets, function(x) all(y %in% x))), ncol = length(TermAtoms))
          rownames(offspring[[term]]) <- names(offspringSets)
          sigOffspring[[term]] <- logical(length(offspringSets))
          names(sigOffspring[[term]]) <- names(offspringSets)
        } else {                                                                # An empty version for when the focus term is an end node
          sigUnions[[term]] <- logical(0)
          atoms[[term]] <- list()
          unions[[term]] <- matrix(,0,0)
          offspring[[term]] <- matrix(,0,0)
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
          basePattern <- offspring[[term]][propagate[[i]],] # The propagated GO term as unions of atoms
          ancestorPatterns <- fillPattern(basePattern)      # All superset patterns as unions of atoms
          ancestorPatterns <- ancestorPatterns[!intersectPatterns(ancestorPatterns, pPatterns),,drop=FALSE]    # Remove duplicates
          pPatterns <- rbind(pPatterns, ancestorPatterns)
        }
        matchedPatterns1 <- intersectPatterns(pPatterns, unions[[term]])  # These new patterns were already testable
        matchedPatterns2 <- intersectPatterns(unions[[term]], pPatterns)  # These already testable patterns are now called significant
        newPatterns <- pPatterns[!matchedPatterns1,, drop=FALSE]          # These are new patterns
        if (verbose) cat("\tSignificance of", nrow(newPatterns), "genesets from", propagate, " propagated to", term, "\n")
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
      offspringPatterns <- offspring[[term]][!sigOffspring[[term]],,drop=FALSE]     # Patterns of not already significant offspring
      matched <- intersectPatterns(offspringPatterns, sigPatterns)                  # Any new offspring terms among the significant patterns?
      newSigOffspring <- rownames(offspringPatterns)[matched]
      if (verbose && (length(setdiff(newSigOffspring, propagate))>0)) {
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
        if (verbose) cat("\ttests:", nrow(newpatterns), "in", term, "\n")
        newsets <- lapply(as.list(1:nrow(newpatterns)), function(i) { # Assemble gene sets from the atoms
          unlist(atoms[[term]][newpatterns[i,]])
        })
        newpvalues <- .globaltest(smallgt, newsets)[,4]
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
    affectedGOforefathers <- unique(unlist(GO@ancestors[newSigGO])) # Find GO terms above the focus level significant through propagation
    affectedGOforefathers <- intersect(affectedGOforefathers, forefathers)
    newlyAffectedGOforefathers <- intersect(affectedGOforefathers, names(adjustedP[is.na(adjustedP)]))
    adjustedP[newlyAffectedGOforefathers] <- alpha
    if (verbose && (length(newlyAffectedGOforefathers) > 0))
      cat("Significant through upward propagation:", newlyAffectedGOforefathers, "\n")
      
    # For all significant GO terms, see if they are also present in other subtrees  
    affectedGOwithParents <- lapply(as.list(newSigGO), function(nsg) { 
      present1 <- sapply(offspring, function(off) nsg %in% rownames(off))
      present2 <- focus %in% nsg
      focus[present1 | present2]
    })      # Creates a list of new significant GO terms, listing the subtrees they appear in
    names(affectedGOwithParents) <- newSigGO
    indirectAffected <- turnListAround(affectedGOwithParents)  # Reverses the list to a list of subtrees, listing the significant GO offspring terms
    
    # Output progress information
    if (!verbose) {
      cat(paste(rep("\b", 48), collapse=""))
      cat("Alpha = ", format(alpha, digits=3, scientific=T, width= 10), ". Significant GO terms: ", format(sum(!is.na(adjustedP)), width=5), ".", sep="")
      flush.console()
    }
    
    # If nothing happened, increase alpha
    if (!change) {                                          
      allPs <- sort(c(pFocus, unlist(pUnions)))
      allPs <- allPs[allPs * holm > alpha]
      alpha <- allPs[1] * holm
      names(alpha) <- NULL
      if (verbose) cat("Alpha =", format(alpha, digits=3, scientific=T, width= 10), "\n")
    } else {
      if (verbose && (holm > sum(!emptyFocus))) cat("Holm's factor:", sum(!emptyFocus), "\n")
      holm <- sum(!emptyFocus)
    }
    ready <- (holm == 0) || (alpha > maxalpha) || (sum(!is.na(adjustedP)) >= stopafter)
  } 
  if (!verbose) cat("\n")
  adjustedP[!is.na(adjustedP)]
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
# Equivalent of %in% but for matrices
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
turnListAround <- function(alist) {
  allcontent <- unique(unlist(alist))
  out <- lapply(as.list(allcontent), function(content) {
    inContainer <- sapply(alist, function(container) content %in% container)
    names(alist[inContainer])
  })
  names(out) <- allcontent
  out
}
