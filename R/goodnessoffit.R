gtPS <- function(response, null, data, model, ..., 
		          covs, bdeg = 3, nint= 10, pord = 2,
				  interact = FALSE, by = NULL, robust = FALSE, termlabels = FALSE) {
  
  # store the call               
  call <- match.call()                	   
                  	            	
  # missing data
  if (missing(data)) data<-NULL
  if (is.matrix(data)) data <- as.data.frame(data)
  
  # evaluate response, which may be one of the colnames of data
  response <- eval(call$response, data, parent.frame())

  # response model
  if (is(response,"lm") | is(response,"glm") | is(response,"coxph")) {
    data <- eval(response$call$data)
     if (is.null(data)) stop("argument \"data\" in the response model is missing") else data<-as.data.frame(data)
    response<-as.formula(response$call$formula)
  }	

  # settle formula null and vector response if response is a formula
  if (missing(null))
    if (is(response, "formula")) 
      null <- response
    else
      stop("argument \"null\" is missing, with no default") 
      
  if (is(response, "formula")) {
    name.response <-  as.character(eval(response)[[2]])
    response <- eval(attr(terms(response, data=data), "variables"), data, environment(response))[[attr(terms(response, data=data), "response")]]
  } else {
    name.response <- deparse(call$response)
  }

  # covs default
  if (missing(covs)) {
  	if (is(null,"formula"))  {
      if (termlabels || (null[[length(null)]]==".")) covs<-attr(terms(null, data=data), "term.labels") 
      else {
        covs<-all.vars(null, functions = FALSE, unique=TRUE)
	    covs<-covs[!(covs%in%name.response)] 
	  }
	  if (robust) X0<-model.matrix(null, data) 
  	}
	if (is.data.frame(null) || is.vector(null) || is.matrix(null)) {
	  null <- data.frame(null)
      intercept <- apply(null,2,function(x){all(x==1)})
        if (all(intercept == FALSE)) covs <- colnames(null)  
        else covs <- colnames(null)[-intercept]
      if (robust) X0<-null
    }
  }

  # build the matrix of selected covs
  if (is(null,"formula")) {
    if (is.null(data)) datacovs<-cbind(model.frame(null,data=data),model.matrix(null,data=data))[,covs,drop=F]
    else datacovs<-cbind(data,model.frame(null,data=data),model.matrix(null,data=data))[,covs,drop=F]
  }
  else datacovs<-null[,covs]
  datacovs<-as.data.frame(as.matrix(datacovs))
  covs<-colnames(datacovs)
  p<-ncol(datacovs)
 
  # interacting factors
  if (!is.null(by)) {
    databy<-data[,by]
    if (!is.null(ncol(databy))) nlev<-lapply(lapply(databy,factor),nlevels) 
    else nlev<-nlevels(factor(databy)) 
  }
  else { databy<-NULL; nlev<-0 }
  
  # print smooth terms
  if (interact) cat("  smooth term:", paste("s(",paste(covs,collapse=","),")", sep=""), "\n") else
  if (!is.null(by)) cat("  smooth terms:", paste(paste("s(",covs,")", sep=""),":",by, sep=""), "\n") 
  else cat("  smooth terms:", paste("s(",covs,")", sep=""), "\n")

  # bdeg|nint|pord is a list with at least two elements
  l <- max(length(bdeg)*is.list(bdeg),length(nint)*is.list(nint),length(pord)*is.list(pord))
  if ((is.list(bdeg)|is.list(nint)|is.list(pord)) && l>1) {
    lnames<-c(names(nint),names(pord),names(bdeg))[1:l]
    lvector<-vector("list",l); names(lvector)<-lnames
	if (!is.list(bdeg))	bdeg<-lapply(lvector,function(x){x<-bdeg})
	if (!is.list(nint))	nint<-lapply(lvector,function(x){x<-nint})
	if (!is.list(pord))	pord<-lapply(lvector,function(x){x<-pord})
	
	#check list
	stopifnot(all(sort(names(bdeg))==sort(names(nint))) && all(sort(names(nint))==sort(names(pord)))) 
	bdeg<-lapply(bdeg,function(x){ if (length(x)==1) rep(x,p) else x})
	nint<-lapply(nint,function(x){ if (length(x)==1) rep(x,p) else x})
	pord<-lapply(pord,function(x){ if (length(x)==1) rep(x,p) else x})
	if (interact) lvector<-lapply(lnames,function(x){prod(unlist(nint[x])+unlist(bdeg[x])) - prod(unlist(pord[x])) }) else 
	lvector<-lapply(lnames,function(x){sum( unlist(nint[x])+unlist(bdeg[x])-unlist(pord[x]) )*max(1,sum(unlist(nlev)))} ) 
	names(lvector)<-lnames
	cvector<-c(0,cumsum(lvector))
	crange<-lapply(1:l,function(x){(cvector[x]+1):(cvector[x+1])}) 
	names(crange)<-lnames
	
	# multiple results
	if (!robust) {
      bigZ<-do.call(cbind,lapply(lapply(lnames, function(x){
      	bbases(datacovs, nint=unlist(nint[x]), bdeg=unlist(bdeg[x]),pord=unlist(pord[x]), interact=interact, by=databy)
      	}),unlist))
	  colnames(bigZ)<-1:ncol(bigZ)
	  res<-gt(response=response, alternative=bigZ, null=null, data=data, ..., subsets=crange )
	} 
	
	# robust result
	else {	
	  H<-X0 %*% solve(crossprod(X0), t(X0))
	  bigZ<-do.call(cbind, lapply( lapply( lapply(lnames, function(x){ 
	  	(diag(nrow(H))-H) %*% bbases(datacovs, nint=unlist(nint[x]), bdeg=unlist(bdeg[x]), pord=unlist(pord[x]), interact=interact, by=databy)
	  	}), function(x){
	  		x*(1/sqrt(sum(diag(t(x) %*% x))))
	  		}), unlist))
      colnames(bigZ)<-1:ncol(bigZ)
      res<-gt(response=response, alternative=bigZ, null=null, data=data, ... )
    }
  }	
  
  # standard result
  else {
	if (length(bdeg)==1) bdeg<-rep(bdeg,p) 
	if (length(nint)==1) nint<-rep(nint,p) 
	if (length(pord)==1) pord<-rep(pord,p) 
	stopifnot(length(c(bdeg,nint,pord))==3*p)
	Z<-bbases(datacovs, bdeg, nint, pord, interact=interact, by=databy)
    res<-gt(response=response, alternative=Z, null=null, data=data, ... )
  }

  return(res) 
}

###################################################
# External functions for gtPS
###################################################

  ### Construct difference matrix of order pord
ndiff <- function(n, pord = 2) {
  if (pord == 0) D <-diag(n)
  else {
  if (pord == 1) { D <- diff(diag(n)) }
  else { D <- diff(ndiff(n, pord - 1)) }
  }
  return(D) 
}
  
  
  ### Truncated p-th power function
tpower <- function(x, t, p)
  (x - t) ^ p * (x > t)


  ### Construct B-spline basis
bbase <- function(x, bdeg, nint){ 
  xl <- min(x) - 0.001*(max(x)-min(x)) 
  xr <- max(x) + 0.001*(max(x)-min(x)) 
  dx <- (xr - xl) / nint
  knots <- seq(xl - bdeg * dx, xr + bdeg * dx, by = dx)
  P <- outer(x, knots, tpower, bdeg)
  n <- dim(P)[2]
  D <- ndiff(n, bdeg + 1) / (gamma(bdeg + 1) * dx ^ bdeg)
  B <- (-1) ^ (bdeg + 1) * P %*% t(D)
  B 
}


  ### Construct tensor product of two or more B-spline bases AND Kronecker sum of the penalty matrices
btp <- function(datacovs, bdeg, nint, pord){
	p <- ncol(datacovs)
	stopifnot( p > 1 )
  	qq <- nint + bdeg
  	stopifnot( prod(qq) < 5000 )	
	bigB<-do.call(cbind, lapply(lapply(1:p,function(x){
		bbase(datacovs[,x], bdeg[x], nint[x])
		}), unlist))
	tpB<-bigB[,( cumsum(qq)[1]-qq[1]+1 ):( cumsum(qq)[1] )]
	ksP<-t(ndiff(ncol(tpB),pord[1])) %*% ndiff(ncol(tpB),pord[1])
	for (j in 2:p) {
      crange<-( cumsum(qq)[j]-qq[j]+1 ):(cumsum(qq)[j])
      qj<-length(crange)
      tpB<-kronecker(tpB, t(rep(1,qj)) ) * kronecker( t(rep(1,ncol(tpB))), bigB[,crange] )
      Pj<-t(ndiff(qj,pord[j])) %*% ndiff(qj,pord[j])
      ksP<-kronecker(diag(qj),ksP) + kronecker( Pj ,diag(ncol(ksP)) )
    }
    return( list(tpB=tpB, ksP=ksP) )
}


  ### get reparameterized design matrices
repdes <-function (B, pord, K=NULL, tol = 1e-10) {
  if (is.null(K)) K = t(ndiff(ncol(B), pord))%*%ndiff(ncol(B), pord)
  ek <- eigen(K, symmetric = TRUE)
  nullvals <- ek$values < tol
  if (any(nullvals)) {
    nullK <- ek$vectors[, nullvals]
    N <- B %*% nullK
    L <- t(sqrt(ek$values[!nullvals]) * t(ek$vectors[, !nullvals]))
    Z <- B %*% L %*% solve(t(L) %*% L)
  }
  else {
    N <- matrix(0, nr = nrow(B), nc = 0)
    Z <- B
  }
  return(list(Z = Z, N = N))
}


  ### Stack reparameterized B-spline bases or get reparameterized tensor product of two or more B-spline bases
bbases <-function(datacovs, bdeg, nint, pord, by, interact) {
  p<-ncol(datacovs)
  
  # interact FALSE
  if (!interact){	
  	qq <- nint + bdeg - pord
    Z<-do.call(cbind, lapply(lapply(1:p,function(x){
    	bbase(datacovs[,x],bdeg[x],nint[x]) %*% t(ndiff(nint[x]+bdeg[x],pord[x])) %*% solve(crossprod(t( ndiff(nint[x]+bdeg[x],pord[x])))) 
    	}),unlist))
    colnames(Z)<-unlist(sapply(1:p,function(x){paste(colnames(datacovs)[x],rep(".",nint[x]+bdeg[x]-pord[x]), sep="") }))
    
	  # by TRUE
	  if (!is.null(by)) {
        nfactors<-max(1,ncol(by))
        # one factor
	    if (nfactors == 1) {
	      f<-factor(by)
		  nlev<-nlevels(f)
		  ZbyF<-do.call(cbind,lapply(lapply(1:p, function(x){
		  	model.matrix(~-1+f:Z[,unlist( (cumsum(qq)[x]-qq[x]+1 ):( cumsum(qq)[x]) )])[,]
		  	}),unlist))
		  ZbyF<-ZbyF[,unlist(lapply(1:nlev,function(x){seq(x,sum(qq)*nlev,nlev)}))]
		  colnames(ZbyF)<-unlist(lapply(1:p,function(x){rep(paste(levels(f),":","s(",colnames(datacovs)[x],")", sep=""),
		  each=(nint[x]+bdeg[x]-pord[x]))}))
		}
		# multiple factors
        else{
		  fs<-lapply(by,factor)
		  nlevs<-unlist(lapply(fs,nlevels))
		  ZbyF<-do.call(cbind,lapply(lapply(1:nfactors, function(x){ 
		  	lapply(1:p, function(xx){ model.matrix(~-1+unlist(fs[x]):Z[,(cumsum(qq)[xx]-qq[xx]+1):(cumsum(qq)[xx])][,])}) 
		  	}),unlist))
	      colnames(ZbyF)<-1:ncol(ZbyF)
	    }
	 Z<-ZbyF   
	 }	
  }
  
  # interact TRUE
  else{
  BK<-btp(datacovs, bdeg=bdeg, nint=nint, pord=pord)
  Z <- repdes(B=BK$tpB, K=BK$ksP)$Z
  colnames(Z)<-1:ncol(Z)
  }
  return(Z)
}


  ### get reparameterized Z 
getPS <-function(data, covs, bdeg=3, nint=10, pord=2, by=NULL, interact=FALSE) {
  datacovs<-data.frame(data[,covs])
  names(datacovs)<-covs
  if (!is.null(dim(datacovs))) {
  	p<-ncol(datacovs) 
    bdeg<-sapply(bdeg,function(x){ if (length(x)==1) rep(x,p) else x})
	nint<-sapply(nint,function(x){ if (length(x)==1) rep(x,p) else x})
	pord<-sapply(pord,function(x){ if (length(x)==1) rep(x,p) else x})
    }
  B<-bbases(datacovs=datacovs, bdeg=bdeg, nint=nint, pord=pord, by=by, interact=interact)
  return(B)
}



gtI <- function(response, null, data, model, ..., covs, iorder=2) {
  
  # match the call	
  call <- match.call()   
  
  # formula response
  formula.response <- (length(call$response) > 1) && (call$response[[1]] == "~")
  if (formula.response) null<- response
  
  # data missing
  if (!missing(data)) data<-as.data.frame(data)
  else data<-model.frame(null) 
  
  # covs default
  if (missing(covs)) {
    if (is.data.frame(null) || is.vector(null) || is.matrix(null)) {
      names.null<-names(as.data.frame(null))
      intercept<-apply(null,2,function(x){all(x==1)})
      if (all(intercept == FALSE)) {
        names.null<-names(as.data.frame(null))
        subfmla <- as.formula( paste("~ 0 +", paste(names.null, collapse= "+")) )
      }
      else {	
        names.null<-names(as.data.frame(null))[-intercept]
        subfmla <- as.formula( paste("~ ", paste(names.null, collapse= "+")) )
      }
      data<-null
      null<-subfmla
    }
    else {
      subfmla<- null 
    }
  }
  
  # specified covs
  else {
    # null matrix
    if (is.data.frame(null) || is.vector(null) || is.matrix(null)) {
      names.null<-names(as.data.frame(null))
      intercept<-apply(null,2,function(x){all(x==1)})
      if (all(intercept == FALSE)) {
        names.null<-names(as.data.frame(null))
        subfmla <- as.formula( paste("~ 0 +", paste(names.null, collapse= "+")) )
      } else {	
        names.null<-names(as.data.frame(null))[-intercept]
        subfmla <- as.formula( paste("~ ", paste(names.null, collapse= "+")) )
      }
      data<-null
      null<-subfmla
    }
    # null formula
    else names.null<-attr(terms(null),"term.labels")
    # check if specified covs are in names.null
    if (!all(covs %in% names.null)) stop("argument \"covs\" does not agree with null terms")
    subfmla<- as.formula( paste("~ ", paste(covs, collapse= "+")) ) 
  }
  
  # construct alternative formula
  interact <-as.formula(paste("~(", paste(attr( terms(subfmla) , "term.labels"), collapse= "+") , ")^",iorder) )
  
  # remove terms from alternative that are also in null
  dup <- attr(terms(interact, data=data), "term.labels") %in% attr(terms(null, data=data), "term.labels")
  if (all(dup)) stop("all covariates in alternative also in null")  
  if (any(dup)) interact <- formula(terms(interact,data=data)[!dup])
  
  # result
  res<-gt(response, alternative=interact, null, data, ...)
  return(res)
}



gtK <- function(response, null, data, model, ..., covs, quant = .25, 
                 metric = c("euclidean", "pearson"), 
                 kernel=c("uniform", "exponential", "triangular","neighbours","Gauss"),
                 robust = FALSE, scale=TRUE, termlabels = FALSE) {

  # store the call               
  call <- match.call()                	   
                  	            	
  # missing data
  if (missing(data)) data<-NULL
  if (is.matrix(data)) data <- as.data.frame(data)
  
  # evaluate response, which may be one of the colnames of data
  response <- eval(call$response, data, parent.frame())

  # response model
  if (is(response,"lm") | is(response,"glm") | is(response,"coxph")) {
    data <- eval(response$call$data)
     if (is.null(data)) stop("argument \"data\" in the response model is missing") else data<-as.data.frame(data)
    response<-as.formula(response$call$formula)
  }	

  # settle half formula null and vector response if response is a formula
  if (missing(null))
    if (is(response, "formula"))
      null <- response
    else
      stop("argument \"null\" is missing, with no default") 

  if (is(response, "formula")) {
    name.response <-  as.character(eval(response)[[2]])
    response <- eval(attr(terms(response, data=data), "variables"), data, environment(response))[[attr(terms(response, data=data), "response")]]
  } else {
    name.response <- deparse(call$response)
  }

  # covs default
  if (missing(covs)) {
  	if (is(null,"formula"))  {
	  if (termlabels || (null[[length(null)]]==".") ) covs<-attr(terms(null, data=data), "term.labels") 
      else {
        covs<-all.vars(null, functions = FALSE, unique=TRUE)
	    covs<-covs[!(covs%in%name.response)] 
      }
	  if (robust) X0<-model.matrix(null, data) 
  	}
	if (is.data.frame(null) || is.vector(null) || is.matrix(null)) {
	  null <- data.frame(null)
      intercept <- apply(null,2,function(x){all(x==1)})
        if (all(intercept == FALSE)) covs <- colnames(null)  
        else covs <- colnames(null)[-intercept]
      if (robust) X0<-null
    }
  }

  # build the matrix of selected covs
  if (is(null,"formula")) {
    if (is.null(data)) datacovs<-cbind(model.frame(null,data=data),model.matrix(null,data=data))[,covs,drop=F]
    else datacovs<-cbind(data,model.frame(null,data=data),model.matrix(null,data=data))[,covs,drop=F]
  }
  else datacovs<-null[,covs]
  datacovs<-as.data.frame(as.matrix(datacovs))
  covs<-colnames(datacovs)
  p<-ncol(datacovs)
  n<-nrow(datacovs)

  # print smooth terms 
  cat("  smooth term:", paste("s(",paste(covs,collapse=","),")", sep=""), "\n") 

  # multiple quantiles 
  l<-length(quant)
  if (l>1) {	
    # robust FALSE
    if (!robust) {
      bigU<-do.call(cbind,lapply(lapply(1:l, function(x){
      	getK(datacovs, covs, quant[x], metric, kernel, scale) 
      	}),unlist))
	  colnames(bigU)<-1:(l*n)
      lvector<-lapply(1:l,function(x){(((x-1)*n)+1):(x*n)})
      names(lvector)<-paste("quant",quant)
	  res<-gt(response=response, alternative=bigU, null=null, data=data, ... ,subsets=lvector)
	}
	# robust TRUE
	else {	
	  H<-X0 %*% solve(crossprod(X0), t(X0)) 
	  bigU<- do.call(cbind,lapply(lapply( lapply(1:l, function(x){ 
	  	(diag(n)-H) %*% getK(datacovs, covs, quant[x], metric, kernel, scale)
	  	} ), function(x){ x * (1/sqrt(sum(diag(t(x) %*% x)))) }) ,unlist))
	  colnames(bigU)<-1:(l*n)
	  res<-gt(response=response, alternative=bigU, null=null, data=data, ... )
	}
  }
  
  # one quantile	
  else {	
  U<-getK(data=datacovs, covs=covs, quant=quant, metric=metric, kernel=kernel, scale=scale) 		
  res<-gt(response=response, alternative=U, null=null, data=data, ... )
  }
  return(res)
}

###################################################
# External functions for gtK
###################################################

getK <- function(data, covs, quant=.25, 
			    metric = c("euclidean", "pearson"), 
				kernel=c("uniform", "exponential", "triangular","neighbours","Gauss"), 
				scale=TRUE ) {
  if (!is.data.frame(data)) stop("data should be a data.frame")
  if (missing(covs)) covs <- names(data)
  if (length(metric)>1) metric <- "euclidean"
  if (!all(covs %in% names(data))) stop("incorrect covariate names")
  metric <- match.arg(metric)
  kernel <- match.arg(kernel)
  
  # metrics (to add: categorical, mixed)
  if ( metric == "euclidean" ) {
  	if (scale) { DD <- as.matrix(dist(scale(data[,covs]))) } 
  	else {  DD <- as.matrix(dist(data[,covs]))  } 
  }
  if ( metric == "pearson" ) {
	DD <- as.matrix(as.dist(1-abs(cor(t(data[,covs])))))
  }
  
  # kernels (to add: Epanechnikov, Biweight, Gauss)
  if ( kernel == "uniform" ) {
    dDD <- DD[diag(nrow(DD)) == 0]
    h <- quantile(dDD, quant)
    U <- apply(DD, c(1,2), function(x) { as.numeric(abs(x) <= h) } )
  }
  if ( kernel == "exponential" ) {
  	dDD <- DD[diag(nrow(DD)) == 0]
    h <- quantile(dDD, quant)
    U <- apply(DD, c(1,2), function(x) { as.numeric(exp((-x)/h)) } )
  }
    if ( kernel == "triangular" ) {
  	dDD <- DD[diag(nrow(DD)) == 0]
    h <- quantile(dDD, quant)
    U <- apply(DD, c(1,2), function(x) { as.numeric( (h-abs(x))*(abs(x)<=h) ) } )
  }
  if ( kernel == "neighbours" ) {
    U <- apply(DD, 1, function(x) { as.numeric(x <= quantile(x, quant)) } )
  }
  if ( kernel == "Gauss" ) {
  	dDD <- DD[diag(nrow(DD)) == 0]
    h <- quantile(dDD, quant)
    U <- apply(DD, c(1,2), function(x) { as.numeric( dnorm(x/h) ) } )
  }
  U <- sweep(U, 2, colSums(U), "/")
  U
}


