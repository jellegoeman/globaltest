inheritance <- function(test, sets, weights, ancestors, offspring,  stop = 1, Shaffer, homogeneous=FALSE, trace) {

  # input checking 1: sets
  if (missing(sets) && is(test, "gt.object") && !is.null(test@subsets)) 
    sets <- test@subsets
  if (missing(sets)) stop("argument \"sets\" is missing with no default")
  if (is.character(sets))
    if (is(test, "gt.object") && !is.null(test@subsets)) {
      test <- test[sets]
      sets <- test@subsets
    }
    if(is(sets,"hclust")) sets=as.dendrogram(sets)
    if(is(sets,"dendrogram"))       sets=dendro2sets(sets)
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
        
  broadstructure=do.broadstructure(sets,ancestors,offspring)
  rm(ancestors,offspring)

  if (missing(weights)) weights <- NULL               # prevents confusion with "weights" method

  if (is(test,"gt.object")) {
    if(missing(Shaffer)) {
      if( length( setdiff(unlist(broadstructure$sets[broadstructure$Shaffernode]),unlist(broadstructure$sets[broadstructure$leaflist])) )==0) Shaffer=TRUE
      else  Shaffer=FALSE }
    else 
      if(Shaffer) if( length( setdiff(unlist(broadstructure$sets[broadstructure$Shaffernode]),unlist(broadstructure$sets[broadstructure$leaflist])) )>0) 
                    stop("There are children of nodes with only leaves as offspring that are not included in the tree-structure. Shaffer can not be applied!")
   }  
  labels= unique(unlist(sets))
  if (missing(trace))    trace <- gt.options()$trace
  if (stop <= 1) {
    maxalpha <- stop
    stopafter <- Inf
  } else {
    maxalpha <- 1
    stopafter <- stop
  }

  # Preprocess gt.object input and calculate focus level p-values
  K <- length(sets)
  test.id=1:K
  nam=names(sets)
  sets=lapply(test.id, function(k) sort(sets[[k]]))
  names(sets)=nam
  rawp=rep(NA,length(sets))
  names(rawp)=nam
  digitsK <- trunc(log10(K))+1
  get.weights<-function() {
    if(is.null(test@subsets)) 
      weights= weights(test)
    else  
      weights= weights(test[which(sapply(test@subsets,function(x){ identical(sort(x),sort(labels) )})[1])])
  }                     
  if (is(test, "gt.object")) {
    if (is.list(test@weights)) if (length(test@weights)>1) stop("The inheritance method is not applicable with individually weighted sets")
    if(missing(weights))    weights= get.weights()
    else if(is.null(weights)) weights= get.weights()
    
    #coumpute weights for all tested subsets
    if(length(weights)!= length(labels)) { if(!is.null(weights)) stop("Weights does not fit the number of leaves.") }
      #extract weigths for leaves nodes
    rawgt <- test
    test <- function(set) rawgt@functions$test(set)[1]
    if(!is.null(rawgt@subsets)) rawgt@subsets=lapply(min(1,length(rawgt@subsets)):length(rawgt@subsets), function(k) sort(rawgt@subsets[[k]]))
    found=sapply(test.id, function(j,target,orig) {which(sapply(1:length(target),function(i,j) identical(target[[i]],orig[[j]]),j ))[1]},rawgt@subsets,sets)
    rawp[!is.na(found)] <- rawgt@result[found[!is.na(found)],1]
    rownames(rawgt@result)[found[!is.na(found)]]= names(sets)[!is.na(found)]
    test.id=setdiff(test.id,which(!is.na(found)))

    result <- t(sapply(test.id, function(i) {
      if (trace) {
        if (i > 1) cat(rep("\b", digitsK+trunc(log10(i-1))+4), sep="")
        cat(i, " / ", K, sep="")
        flush.console()
      }
      rawgt@functions$test(sets[[i]])
    }))
    if (length(test.id)>0) {
      rawp[test.id]= result[,1]
      rawgt@result <- rbind(rawgt@result[found[!is.na(found)],],result)
    }        
    rownames(rawgt@result)<- c(names(sets)[!is.na(found)], names(sets)[which(is.na(found))])
    colnames(rawgt@result) <- c("p-value", "Statistic", "Expected", "Std.dev", "#Cov")
    rawgt@subsets <- sets
    if (any(is.na(found))) {
      alias <- alias(rawgt)
      rawgt@extra <- NULL
      if (!is.null(alias)) 
        alias(rawgt) <- c(alias[found[!is.na(found)]], rep("", length(rawgt)-length(alias)))
    }
  } else {
    rawgt <- NULL
    rawp <- sapply(1:K, function(i) {
      if (trace) {
        cat(rep("\b", 2*digitsK+3), i, " / ", K, sep="")
        flush.console()
      }
      test(sets[[i]])
    })
    if(missing(weights)) weights=rep(1,length(rawp))
    else if(is.null(weights)) weights=rep(1,length(rawp))
  }
  if (trace) cat(rep("\b", 2*digitsK+25), sep="")
  #extract weigths for leaves nodes
  weights.sets=sapply(sets,function(x) sum(weights[x]))
  #now recall core function:
  res=.inherit(ps=rawp,structure=broadstructure, weights=weights.sets,stop=maxalpha,Shaffer=Shaffer,homogeneous)
  adjustedP=res$adj.p[names(rawgt)]
  if (trace==1) cat("\n")
  if (!is.null(rawgt)) {
    extra <- rawgt@extra
    extra[["inheritance"]][match(names(adjustedP),names(rawgt))] <- adjustedP
    rawgt@extra <- as.data.frame(extra) 
    rownames(rawgt@extra)=names(rawgt)
#    rawgt=rawgt[sort.list(row.names(rawgt@result)),]
    rawgt@structure <- list(ancestors=broadstructure$ancestor,offspring=broadstructure$offspring) 
    return(rawgt)
  } else {
    return(data.frame(raw.p = rawp, inheritance = adjustedP))
  }

}

#############################################################################################################
################################################# FORMAT CONVERTING FUNCTIONS, SOME OTHER USEFUL FUNCTIONS and the core function

##################################################
do.broadstructure=function(sets,ancestor,offspring){
####################
  is.subset<-  function(sups,subs){ (sum(sets[[subs]] %in% sets[[sups]])==length(sets[[subs]]))&(length(sets[[sups]])>length(sets[[subs]])) }
  is.supset<-  function(subs,sups){ (sum(sets[[subs]] %in% sets[[sups]])==length(sets[[subs]]))&(length(sets[[sups]])>length(sets[[subs]])) }
####################

   parent=ancestors2parents(ancestor)
   child=ancestors2parents(offspring)
   leaf.list<-function(structure) setdiff(names(structure$ancestor),names(structure$offspring))
   ll= leaf.list(list(ancestor=ancestor ,offspring=offspring))
   lapply(child, function(leaf) length(intersect(leaf[[1]],ll))==length(leaf[[1]])  )
    Shaffernode=lapply(child, function(leaf){  length(intersect(leaf,ll))==length(leaf) } )
    Shaffernode=names(Shaffernode)[unlist(Shaffernode)]

  list(ancestor=ancestor ,offspring=offspring,sets=sets, parent=parent,child=child, leaflist=ll, Shaffernode=Shaffernode)
}

################################################################################
dendro2sets<-function(hc){
  ####################
  #create set names from the dendrogram
  do.sets<-function(x,parent.labels="",sets){
    for(i in 1:length(x)){
        ans=list(labels(x[[i]]))
        names(ans)=paste(parent.labels,"[",i,sep="")
        if(!is.null(attributes(x[[i]])$label)){
          names(ans)=paste(names(ans),":",attributes(x[[i]])$label,sep="") }
        sets=c(sets,ans)
        if(is.null(attributes(x[[i]])$leaf)) {sets=do.sets(x=x[[i]],parent.labels=names(ans),sets)  }
    }
  return=sets}
    ####################

  sets=list();
  sets[[1]]=labels(hc)
  names(sets)[1]="O1"
  sets=do.sets(hc,names(sets)[1],sets)          
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
#################################################
# Takes sets to ancestors mapping and converts it to a sets of parents (as well offspring to child)
ancestors2parents <- function(ancestors) {
  lapply(ancestors, function(anc) {
    setdiff(anc, unique(unlist(ancestors[anc])))
  })
} #see the original function in focuslevel.R

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
      


######################################################################################################################
## here start the main function
.inherit<-function(ps,structure, weights,stop=1,Shaffer=T,homogeneous=F){
                                                              
  ps=signif(ps,digits =8)
  m=length(structure$leaflist)
  if(length(weights)==1)  weights=rep(weights,m)
  adj.p=rep(1,length(ps))
  names(adj.p)=names(ps)

  alphas=rep(0,length(ps)); names(alphas)=names(ps)
  founders=setdiff(names(ps),names(structure$ancestor))
  who.next=which.min(ps[founders]/weights[founders])
  alpha=signif(ps[founders[who.next]]/weights[founders[who.next]]*sum(weights[founders]),digits=8)
  alphas[founders]=weights[founders]/sum(weights[founders])*alpha
                                                   
  actives=rep(FALSE,length(alphas)); names(actives)=names(ps)
  rejected=actives
  actives[founders]=TRUE
  istep= TRUE
  reject.previous=0
  while(!all(rejected)){  #repeat until max alpha is reached or all hypotheses are already rejected

    #### start che Inheritance procedure with fixed alpha and procede until as much rejections as possible are done!!
    while( istep ){
      #### start Shaffer
      #### compute Shaffer improvement if at least one lead node is tested
      alphas.Shaffer=alphas
      if(Shaffer & (any(actives[structure$leaflist]))){ #if Shaffer=T and not all leaf nodes are already rejected or actives
        activesleaves=structure$leaflist[actives[structure$leaflist]]
        while(length(activesleaves)){
          sibs=structure$child[[structure$parent[[activesleaves[1]]]]]
		      if (all(sibs %in% activesleaves)) { 
            alphas.Shaffer[sibs]=alphas[sibs]*sum(weights[sibs])/(sum(weights[sibs])-min(weights[sibs])) 
          }
		      activesleaves=setdiff(activesleaves,sibs)
        }
      }
      #### end Shaffer
      
      #### do the following for all actives nodes
      for(i.divide in names(which(actives)) ){
#	  print(c(i.divide,ps[i.divide],alphas.Shaffer[i.divide]))
        if( ps[i.divide]<=alphas.Shaffer[i.divide]){ #if rejected
          rejected[i.divide]=TRUE
          actives[i.divide]=FALSE
          to.divide=alphas[i.divide]
          if(!(i.divide%in%structure$leaflist)){ #if not leaf node   (and rejected)
            ids=structure$child[[i.divide]]
            actives[ids]=TRUE
            #redistribute alpha
            alphas[ ids ]=to.divide*weights[ids]/sum(weights[ids])
            alphas[i.divide]=0
          }
          else{ #if a lead node    (and rejected)   
           if(homogeneous){  ####ALTERNATIVE WAY: redistribute EQUALLY the alphas left free from rejected alphas
             alphas[i.divide]=0
             alphas=alphas*alpha/(alpha-to.divide)
           }else{ #STANDARD inheritance method
           ##### redistribute the free alpha to the heirs (the "inheriters")
            search.heirs=TRUE
            target=i.divide
            while(search.heirs){   
               if(is.null(structure$parent[[target]])){ search.heirs=FALSE}
               else{
                   heirs=setdiff(names(which(actives[structure$offspring[[structure$parent[[target]]]]])),target)
                   if(length(heirs)){
                     search.heirs=FALSE
                     actives[heirs]=TRUE        
                     if (sum(weights[heirs]) > 0)
                       alphas[heirs]=alphas[heirs]+to.divide*weights[heirs]/sum(weights[heirs])
                     alphas[i.divide]=0
                   }
                   else{ target=structure$parent[[target]]}
                }
             }            
           }##<-- END    stop research  
         }  ## end of "if a lead node    (and rejected)"
        } #end if rejected
      } #end for loop
      istep= sum(rejected)>reject.previous
      reject.previous=sum(rejected)
    }    #end while
   #assign  adjusted p-values 
   adj.p[rejected]=sapply(adj.p[rejected],min,alpha,na.rm=TRUE)
   ##increase alpha
   if(!all(rejected)){
      test.next=names(ps)[actives]
      who.next=test.next[which.min(ps[test.next]/weights[test.next])]
      alphas[test.next]=ps[who.next]/alphas[who.next]  *alphas[test.next]
      alpha=sum(alphas[test.next])
      if (is.na(alpha)) 
        alpha <- max(1, stop)
      istep=TRUE
      if(alpha>=stop) rejected[]=TRUE
    }
  }  #end while min p-value is reached

  return(list(adj.p=unlist(adj.p),p.value=ps,structure=structure))
}    #end function