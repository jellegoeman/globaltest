inheritance <- function(test, hc, w, stop = 1, Shaffer, homogeneous=FALSE, trace) {

  if (missing(hc)) {
    if (is(test,"gt.object")) {
      if(!is.null(slot(test, "structure")))
        {structure=test@structure          
         sets=unique(unlist(structure))
         if(is.null(structure$ancestor)|is.null(structure$offspring))  {
            broadstructure=do.structure(sets)
            structure=broadstructure[c("ancestor","offspring")]      }
         else{
            parent=ancestors2parents(structure$ancestors)
            broadstructure=c(structure,list(sets, parent=ancestors2parents(structure$ancestors))) }}
      else stop("argument \"hc\" is missing with no default")
    }
  }
  else{
    if(is(hc,"hclust")) hc=as.dendrogram(hc)
    if(is(hc,"dendrogram")) {
      broadstructure=dendro2structure(hc)
      sets=broadstructure[["sets"]]
      structure=broadstructure[c("ancestor","offspring")] }
    else {if(is(hc,"list"))  {sets=hc; broadstructure=do.structure(sets);   structure=broadstructure[1:2]}
          else  stop("Format of argument \"hc\" is not \"dendrogram\" nor \"hclust\" nor a list of sets")}
  }

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
  if (trace) cat("Getting raw p-values: ")
  K <- length(sets)
  test.id=1:K
  nam=names(sets)
  sets=lapply(test.id, function(k) sort(sets[[k]]))
  names(sets)=nam
  rawp=rep(NA,length(sets))
  names(rawp)=nam
  digitsK <- trunc(log10(K))+1

  if (is(test, "gt.object")) {
    if (is.list(test@weights)) if (length(test@weights)>1) stop("The inheritance method is not applicable with individually weighted sets")
         if(missing(w)) w=test@weights
         else if(is.null(w)) w=test@weights
    #coumpute weights for all tested subsets
    if(length(w)!= length(labels)) { if(!is.null(w)) print("Weights does not fit the number of univariate hypotheses. Uniform weights are assumed") ; w=rep(1,length(labels)) ; names(w)=labels; }
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
     rawp[test.id]= result[,1]

     rawgt@result <- rbind(rawgt@result[found[!is.na(found)],],result)
     rownames(rawgt@result)<- c(names(sets)[!is.na(found)], names(sets)[which(is.na(found))])
     colnames(rawgt@result) <- c("p-value", "Statistic", "Expected", "Std.dev", "#Cov")
     rawgt@subsets <- sets
     rawgt@structure=structure
  } else {
    rawgt <- NULL
    rawp <- sapply(1:K, function(i) {
      if (trace) {
        cat(rep("\b", 2*digitsK+3), i, " / ", K, sep="")
        flush.console()
      }
      test(sets[[i]])
    })
    if(missing(w)) w=rep(1,length(rawp))
    else if(is.null(w)) w=rep(1,length(rawp))
  }
  if (trace) cat(rep("\b", 2*digitsK+25), sep="")
  #extract weigths for leaves nodes
  weights=sapply(sets,function(x) sum(w[x]))
  #now recall core function:
  res=.inherit(ps=rawp,structure=broadstructure, w=weights,stop=maxalpha,Shaffer=Shaffer,homogeneous)
  adjustedP=res$adj.p[names(rawgt)]
  if (trace==1) cat("\n")
  if (!is.null(rawgt)) {
    extra <- rawgt@extra
    extra[["inheritance"]][match(names(adjustedP),names(rawgt))] <- adjustedP
    rawgt@extra <- as.data.frame(extra)
    rownames(rawgt@extra)=names(rawgt)
    #rawgt=rawgt[sort.list(row.names(rawgt@result)),] 
    return(rawgt)
  } else {
    return(data.frame(raw.p = rawp, inheritance = adjustedP))
  }

}

#############################################################################################################
################################################# FORMAT CONVERTING FUNCTIONS, SOME OTHER USEFUL FUNCTIONS and the core function
do.structure=function(sets){
####################
  is.subset<-  function(sups,subs){ (sum(sets[[subs]] %in% sets[[sups]])==length(sets[[subs]]))&(length(sets[[sups]])>length(sets[[subs]])) }
  is.supset<-  function(subs,sups){ (sum(sets[[subs]] %in% sets[[sups]])==length(sets[[subs]]))&(length(sets[[sups]])>length(sets[[subs]])) }
####################
  ancestor=list()
  offspring=list()
  for(i in 1:length(sets)){
   temp=names(sets)[sapply(1:length(sets), is.subset,i)]
   if(length(temp)>0){
     ancestor=c(ancestor,list(temp))
     names(ancestor)[length(ancestor)]=names(sets)[i]
     names(sets)[i]
     }
   temp=names(sets)[sapply(1:length(sets), is.supset,i)]
   if(length(temp)>0){
     offspring=c(offspring,list(temp))
     names(offspring)[length(offspring)]= names(sets)[i]
     names(sets)[i]
     }
   }
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
dendro2structure<-function(hc){
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
  ans=do.structure(sets) 
}


#################################################
# returns the parent of element id
get.parent<-function(id){
  temp=sapply(structure$offspring[structure$ancestor[[id]]],length)
  structure$ancestor[[id]][temp == min(c(unlist(temp),Inf))]
  }
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
# Takes sets to ancestors mapping; converts it to sets to parents mapping
ancestors2parents <- function(ancestors) {
  lapply(ancestors, function(anc) {
    setdiff(anc, unique(unlist(ancestors[anc])))
  })
} #see the original function in focuslevel.R

######################################################################################################################
## here start the main function
.inherit<-function(ps,structure=structure, w,stop=1,Shaffer=T,homogeneous=F){
  ps=signif(ps,digits =8)
  m=length(structure$leaflist)
  if(length(w)==1)  w=rep(w,m)
  adj.p=rep(1,length(ps))
  names(adj.p)=names(ps)

  alphas=rep(0,length(ps)); names(alphas)=names(ps)
  founders=setdiff(names(ps),names(structure$ancestor))
  who.next=which.min(ps[founders]/w[founders])
  alpha=signif(ps[founders[who.next]]/w[founders[who.next]]*sum(w[founders]),digits=8)
  alphas[founders]=w[founders]/sum(w[founders])*alpha

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
		       if(all(sibs %in% activesleaves)) { alphas.Shaffer[sibs]=alphas[sibs]*sum(w[sibs])/(sum(w[sibs])-min(w[sibs])) }
		       activesleaves=setdiff(activesleaves,sibs)
          }
      }
      #### end Shaffer
      
      #### do the following for all actives nodes
      for(i.divide in names(which(actives)) ){
        if( ps[i.divide]<=alphas.Shaffer[i.divide]){ #if rejected
          rejected[i.divide]=TRUE
          actives[i.divide]=FALSE
          to.divide=alphas[i.divide]
          if(!(i.divide%in%structure$leaflist)){ #if not leaf node   (and rejected)
            ids=structure$child[[i.divide]]
            actives[ids]=TRUE
            #redistribute alpha
            alphas[ ids ]=to.divide*w[ids]/sum(w[ids])
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
                     alphas[heirs]=alphas[heirs]+to.divide*w[heirs]/sum(w[heirs])
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
      who.next=test.next[which.min(ps[test.next]/w[test.next])]
      alphas[test.next]=ps[who.next]/alphas[who.next]  *alphas[test.next]
      alpha=sum(alphas[test.next])
      istep=TRUE
      if(alpha>stop) rejected[]=TRUE
    }
  }  #end while min p-value is reached

return(list(adj.p=unlist(adj.p),p.value=ps,structure=structure))
}    #end function