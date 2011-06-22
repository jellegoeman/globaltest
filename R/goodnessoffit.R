gtPS <- function(response, null, data, model = c("linear", "logistic", "cox", "poisson", "multinomial"), ..., 
  	          covs, bdeg = 3, nint= 10, pord = 2, interact = FALSE, robust = FALSE, 
				      termlabels = FALSE, returnZ = FALSE) {
  
  # store the call               
  call <- match.call()                	   
                  	            	
  # missing data
  if (missing(data)) data<-NULL
  if (is.matrix(data)) data <- as.data.frame(data)
  
  # evaluate response, which may be one of the colnames of data
  response <- eval(call$response, data, parent.frame())

  # response model
  if (is(response,"lm") || is(response,"glm") || is(response,"coxph")) {
    data <- eval(response$call$data)
    if (is.null(data)) stop("argument \"data\" in the null model is missing") 
    else data<-as.data.frame(data)
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
  	}
	if (is.data.frame(null) || is.vector(null) || is.matrix(null)) {
	  null <- data.frame(null)
    intercept <- apply(null,2,function(x){all(x==1)})
      if (all(intercept == FALSE)) covs <- colnames(null)  
      else covs <- colnames(null)[-intercept]
    }
  }

  # get the model
  if (missing(model)) {
    if(is(response, "Surv")) 
      model <- "cox"
    else if ((is.factor(response) && length(levels(response)) <= 2) || is.logical(response)) 
      model <- "logistic"
    else if (is.factor(response) && length(levels(response)) > 2)
      model <- "multinomial"
  }
  model <- match.arg(tolower(model), c("linear", "logistic", "cox", "poisson", "multinomial"))

  # null design X0
  if (is(null, "formula")) X0<-model.matrix(null, data) else X0<-null
  m <- ncol(X0)
  n <- nrow(X0)
  form <- formula("response ~ 0 + X0")
  if (model=="cox"){
    intercept <- apply(X0,2,function(x){all(x==1)})
    X0<- X0[,-intercept, drop=F]
  }
  
  # make the hat matrix H
  if (model=="linear") 
    H <- X0 %*% solve(crossprod(X0), t(X0))
  else {
    null.fit <- switch(model,
      logistic = glm(form,  family = "binomial", data=data),
      poisson = glm(form, family = "poisson", data=data),
      cox = coxph(form, method = "breslow", data=data),
      multinomial = mlogit(form, data=data)
    )
    if (model=="logistic" || model=="poisson") {
      W <- null.fit$weights
      sqrtW <- sqrt(W)
      X0Whf <- X0 * matrix(sqrtW, n, m)
      X0W <- X0 * matrix(W, n, m)
      H <- X0 %*% solve(crossprod(X0Whf), t(X0W))
     }
    else if (model == "cox") {
      expci <- exp(null.fit$linear.predictors)
      times <- as.vector(null.fit$y)[1:n]
      d <- as.vector(null.fit$y)[(n+1):(2*n)]
      dtimes <- unique(times[d == 1])
      nd <- length(dtimes)
      matrixO <- outer(times, dtimes, "==") * matrix(d, n, nd)
      ties <- any(colSums(matrixO) > 1)
      atrisk <- outer(times, dtimes, ">=")
      hazinc <- as.vector(1 / (expci %*% atrisk))
      matrixP <- outer(expci, hazinc, "*") * atrisk
      matrixPO <- matrixP %*% t(matrixO)
      matrixM <- matrix(d, n, nd, byrow=FALSE) * (!atrisk) - matrixPO %*% (!atrisk)
      matrixW <- -crossprod(t(matrixPO))
      diag(matrixW) <- diag(matrixW) + rowSums(matrixPO)
      X0W <- crossprod(X0, matrixW)
      X0WX0 <- X0W %*% X0
      X0WX0invX0 <- solve(X0WX0, t(X0))
      H <- X0 %*% solve(X0WX0, X0W)
    }
    else if (model == "multinomial") { #temporary
      H <- X0 %*% solve(crossprod(X0), t(X0))
    }
  }
  
  # smooth terms
  if (interact) sterms<-paste(paste("s(",paste(covs,collapse=","),")", sep=""), covs)
  else sterms<-paste("s(",covs,")", sep="")
  
  # build datacovs
  if (is(null,"formula")) {
    if (is.null(data)) datacovs<-cbind(model.frame(null,data=data),model.matrix(null,data=data))[,covs,drop=F]
    else datacovs<-cbind(data,model.frame(null,data=data),model.matrix(null,data=data))[,covs,drop=F]
  }
  else datacovs<-null[,covs]
  datacovs<-as.data.frame(as.matrix(datacovs))
  covs<-colnames(datacovs)
  p<-ncol(datacovs)

  # arguments: list
  l <- max(length(bdeg)*is.list(bdeg),length(nint)*is.list(nint),length(pord)*is.list(pord))
  if ((is.list(bdeg)|is.list(nint)|is.list(pord)) && l>1) {
    lnames<-c(names(nint),names(pord),names(bdeg))[1:l]
    lvector<-vector("list",l); names(lvector)<-lnames
	if (!is.list(bdeg))	bdeg<-lapply(lvector,function(x){x<-bdeg})
	if (!is.list(nint))	nint<-lapply(lvector,function(x){x<-nint})
	if (!is.list(pord))	pord<-lapply(lvector,function(x){x<-pord})
	
	# check the list of covs
	stopifnot(all(sort(names(bdeg))==sort(names(nint))) && all(sort(names(nint))==sort(names(pord)))) 
	bdeg<-lapply(bdeg,function(x){ if (length(x)==1) rep(x,p) else x})
	nint<-lapply(nint,function(x){ if (length(x)==1) rep(x,p) else x})
	pord<-lapply(pord,function(x){ if (length(x)==1) rep(x,p) else x})
	if (interact) lvector<-lapply(lnames,function(x){prod(unlist(nint[x])+unlist(bdeg[x])) - prod(unlist(pord[x])) }) else 
	lvector<-lapply(lnames,function(x){sum( unlist(nint[x])+unlist(bdeg[x])-unlist(pord[x]) ) } ) 
	names(lvector)<-lnames
	cvector<-c(0,cumsum(lvector))
	crange<-lapply(1:l,function(x){(cvector[x]+1):(cvector[x+1])}) 
	names(crange)<-lnames
	
	# built Z matrix
	if (interact){
	Z<-do.call(cbind, lapply(lapply(lapply(lnames, function(x){ 
		(diag(nrow(H))-H) %*%  btensor(datacovs, unlist(bdeg[x]), unlist(nint[x]), unlist(pord[x])) 
		}), function(z){
	  		z*(1/sqrt(sum(diag(crossprod(z)))))
      		}),unlist))
      	colnames(Z)<-1:ncol(Z)
    }
	else {
	Z<-do.call(cbind, lapply(
	lapply(lnames, function(x){
	do.call(cbind, lapply( lapply( lapply(1:p, function(y){
     (diag(nrow(H))-H) %*% bbase(datacovs[,y],unlist(bdeg[x])[y],unlist(nint[x])[y]) %*% 
     t(ndiff(unlist(nint[x])[y]+unlist(bdeg[x])[y],unlist(pord[x])[y])) %*% 
     solve(crossprod(t( ndiff(unlist(nint[x])[y]+unlist(bdeg[x])[y],unlist(pord[x])[y])))) 
     }), function(z){
	 	z*(1/sqrt(sum(diag(crossprod(z)))))
      	}),unlist))
	}), unlist))
	colnames(Z)<-unlist(lapply(lnames, function(y){
		sapply(1:p,function(x){rep(colnames(datacovs)[x],(unlist(nint[y])[x]+unlist(bdeg[y])[x]-unlist(pord[y])[x])) })
	}))
	}
  
  # result
	if (!robust) res<-gt(response=response, alternative=Z, null=null, data=data, model=model,..., subsets=crange )
	else res<-gt(response=response, alternative=Z, null=null, data=data, model=model, ... )
	res@structure <- lapply(lnames,function(x){data.frame(s.term=sterms,bdeg=unlist(bdeg[x]),nint=unlist(nint[x]),pord=unlist(pord[x]))})
  }	
  
  # arguments: vector
  else {
    if (is.list(bdeg)|is.list(nint)|is.list(pord)) {bdeg<-unlist(bdeg);nint<-unlist(nint);pord<-unlist(pord)}
	  if (length(bdeg)==1) bdeg<-rep(bdeg,p) 
	  if (length(nint)==1) nint<-rep(nint,p) 
	  if (length(pord)==1) pord<-rep(pord,p) 
	  stopifnot(length(c(bdeg,nint,pord))==3*p)
	
  # build Z
	  if (interact) {
		  Z<-btensor(datacovs, bdeg, nint, pord)
	    colnames(Z)<-1:ncol(Z)
	  }
	  else {
	    Z<-do.call(cbind, lapply( lapply( lapply(1:p, function(x){
             (diag(nrow(H))-H) %*% bbase(datacovs[,x],bdeg[x],nint[x]) %*% t(ndiff(nint[x]+bdeg[x],pord[x])) %*% solve(crossprod(t( ndiff(nint[x]+bdeg[x],pord[x])))) 
              }),
      	 function(y){
	  		y*(1/sqrt(sum(diag(crossprod(y)))))
      	}),unlist))
      colnames(Z)<-unlist(sapply(1:p,function(x){paste(colnames(datacovs)[x],1:(nint[x]+bdeg[x]-pord[x]), sep="") }))
	  }
  
  # result
  res<-gt(response=response, alternative=Z, null=null, data=data, model=model,... )
  res@structure <- list(data.frame(s.term=sterms,bdeg,nint,pord))    
  }

  if (!returnZ) return(res)
  else return(list(result=res, Z=Z))
}

###################################################
# External functions for gtPS
###################################################

  ### Extract smooth terms from a gt.object object
setGeneric("sterms", function(object, ...) standardGeneric("sterms"))
setMethod("sterms", "gt.object",
  function(object) {
  if (length(object@structure)==1) as.data.frame(object@structure)
  else object@structure
  }
)


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


  ### Construct tensor product of B-spline bases (reparameterized by Kroneker sum of penalties)
btensor <- function(xs, bdeg, nint, pord, returnU=FALSE){
	p <- ncol(xs)
	stopifnot( p > 1 )
  	qq <- nint + bdeg
  	stopifnot( prod(qq) < 5000 )	
	bigB<-do.call(cbind, lapply(lapply(1:p,function(x){
		bbase(xs[,x], bdeg[x], nint[x])
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
    P <- reparamZ(Z=tpB, K=ksP)
    if (!returnU) return(P)
    else  {
    U <- reparamZ(Z=tpB, K=ksP, returnU=TRUE)$U
    return(list(P=P, U=U)) }
}


  ### Reparameterize Z
reparamZ <-function (Z, pord, K=NULL, tol = 1e-10, returnU=FALSE) {
  if (is.null(K)) K = t(ndiff(ncol(Z), pord)) %*% ndiff(ncol(Z), pord)
  ek <- eigen(K, symmetric = TRUE)
  nullvals <- ek$values < tol
  if (any(nullvals)) {
    nullK <- ek$vectors[, nullvals]
    U <- Z %*% nullK
    colnames(U)<-1:ncol(U)
    L <- t(sqrt(ek$values[!nullvals]) * t(ek$vectors[, !nullvals]))
    P <- Z %*% L %*% solve(t(L) %*% L)
    colnames(P)<-1:ncol(P)
  }
  else {
    U <- matrix(0, nr = nrow(Z), nc = 0)
    P <- Z
  }
  if (!returnU) return(P)
  else return(list(P=P, U=U))
}

  ### Standardize Z
reweighZ <- function(Z, null.fit){
  X0 <- model.matrix(null.fit)
  m <- ncol(X0)
  n <- nrow(X0)
  # make the hat matrix H
  if (attr(null.fit,"class")[1]=="lm") 
    H <- X0 %*% solve(crossprod(X0), t(X0))
  else if ( (attr(null.fit,"class")[1]=="glm")  && (null.fit$family[1]=="binomial" || null.fit$family[1]=="poisson") ) {
    W <- null.fit$weights
    sqrtW <- sqrt(W)
    X0Whf <- X0 * matrix(sqrtW, n, m)
    X0W <- X0 * matrix(W, n, m)
    H <- X0 %*% solve(crossprod(X0Whf), t(X0W))
  }
  else if (attr(null.fit,"class")[1]=="coxph") {
    expci <- exp(null.fit$linear.predictors)
    times <- as.vector(null.fit$y)[1:n]
    d <- as.vector(null.fit$y)[(n+1):(2*n)]
    dtimes <- unique(times[d == 1])
    nd <- length(dtimes)
    matrixO <- outer(times, dtimes, "==") * matrix(d, n, nd)
    ties <- any(colSums(matrixO) > 1)
    atrisk <- outer(times, dtimes, ">=")
    hazinc <- as.vector(1 / (expci %*% atrisk))
    matrixP <- outer(expci, hazinc, "*") * atrisk
    matrixPO <- matrixP %*% t(matrixO)
    matrixM <- matrix(d, n, nd, byrow=FALSE) * (!atrisk) - matrixPO %*% (!atrisk)
    matrixW <- -crossprod(t(matrixPO))
    diag(matrixW) <- diag(matrixW) + rowSums(matrixPO)
    X0W <- crossprod(X0, matrixW)
    X0WX0 <- X0W %*% X0
    X0WX0invX0 <- solve(X0WX0, t(X0))
    H <- X0 %*% solve(X0WX0, X0W)
    }    
  Zadj<-(diag(nrow(H))-H) %*% Z
  Zstar <- Zadj*(1/sqrt(sum(diag(crossprod(Zadj)))))
  return(Zstar)
}


gtLI <- function(response, null, data, ..., covs, iorder=2,
                termlabels = FALSE, standardize = FALSE) {
  
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
  	}
	if (is.data.frame(null) || is.vector(null) || is.matrix(null)) {
	  null <- data.frame(null)
    intercept <- apply(null,2,function(x){all(x==1)})
      if (all(intercept == FALSE)) covs <- colnames(null)  
      else covs <- colnames(null)[-intercept]
    }
  }
  
  covsformula<- as.formula( paste("~ ", paste(covs, collapse= "+")) ) 
  
  # construct alternative formula
  interactions <-as.formula(paste("~(", paste(attr( terms(covsformula) , "term.labels"), collapse= "+") , ")^",iorder) )
  
  # remove terms from alternative that are also in null 
  dup <- attr(terms(interactions, data=data), "term.labels") %in% attr(terms(covsformula, data=data), "term.labels")
  if (all(dup)) stop("all covariates in alternative also in null")  
  if (any(dup)) interactions <- formula(terms(interactions,data=data)[!dup])
  if (is(null,"formula"))  {
  dup2 <- attr(terms(interactions, data=data), "term.labels") %in% attr(terms(null, data=data), "term.labels")
  if (all(dup2)) stop("all covariates in the alternative are also in the null")  
  if (any(dup2)) interactions <- formula(terms(interactions,data=data)[!dup])
  }
  
  # result
  res<-gt(response=response, alternative=interactions, null=null, data=data, standardize=standardize, ...)
  res@structure<-list(interactions=attr(terms.formula(interactions),"term.labels"))
  return(res)
}


gtKS <- function(response, null, data, model = c("linear", "logistic", "cox", "poisson", "multinomial"), ..., 
                covs, quant = .25, metric = c("euclidean", "pearson"), 
                kernel=c("uniform", "exponential", "triangular","neighbours","gauss"),
                robust = FALSE, scale=TRUE, termlabels = FALSE, returnZ = FALSE) {

  # store the call               
  call <- match.call()
  metric <- match.arg(metric)
  kernel <- match.arg(kernel)                	   
   	
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
  	}
	if (is.data.frame(null) || is.vector(null) || is.matrix(null)) {
	  null <- data.frame(null)
      intercept <- apply(null,2,function(x){all(x==1)})
        if (all(intercept == FALSE)) covs <- colnames(null)  
        else covs <- colnames(null)[-intercept]
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

 # get the model
  if (missing(model)) {
    if(is(response, "Surv")) 
      model <- "cox"
    else if ((is.factor(response) && length(levels(response)) <= 2) || is.logical(response)) 
      model <- "logistic"
    else if (is.factor(response) && length(levels(response)) > 2)
      model <- "multinomial"
  }
  model <- match.arg(tolower(model), c("linear", "logistic", "cox", "poisson", "multinomial"))

  # null design X0
  if (is(null, "formula")) X0<-model.matrix(null, data) else X0<-null
  m <- ncol(X0)
  n <- nrow(X0)
  form <- formula("response ~ 0 + X0")
    
  # make the hat matrix H
  if (model=="linear") 
    H <- X0 %*% solve(crossprod(X0), t(X0))
  else {
    null.fit <- switch(model,
      logistic = glm(form, family = "binomial", data=data),
      poisson = glm(form, family = "poisson", data=data),
      cox = coxph(form, method = "breslow", data=data),
      multinomial = mlogit(form, data=data)
    )
    if (model=="logistic" || model=="poisson") {
      W <- null.fit$weights
      sqrtW <- sqrt(W)
      X0Whf <- X0 * matrix(sqrtW, n, m)
      X0W <- X0 * matrix(W, n, m)
      H <- X0 %*% solve(crossprod(X0Whf), t(X0W))
     }
    else if (model == "cox") {
      expci <- exp(null.fit$linear.predictors)
      times <- as.vector(null.fit$y)[1:n]
      d <- as.vector(null.fit$y)[(n+1):(2*n)]
      dtimes <- unique(times[d == 1])
      nd <- length(dtimes)
      matrixO <- outer(times, dtimes, "==") * matrix(d, n, nd)
      ties <- any(colSums(matrixO) > 1)
      atrisk <- outer(times, dtimes, ">=")
      hazinc <- as.vector(1 / (expci %*% atrisk))
      matrixP <- outer(expci, hazinc, "*") * atrisk
      matrixPO <- matrixP %*% t(matrixO)
      matrixM <- matrix(d, n, nd, byrow=FALSE) * (!atrisk) - matrixPO %*% (!atrisk)
      matrixW <- -crossprod(t(matrixPO))
      diag(matrixW) <- diag(matrixW) + rowSums(matrixPO)
      X0W <- crossprod(X0, matrixW)
      X0WX0 <- X0W %*% X0
      X0WX0invX0 <- solve(X0WX0, t(X0))
      H <- X0 %*% solve(X0WX0, X0W)
    }
    else if (model == "multinomial") {
      H <- X0 %*% solve(crossprod(X0), t(X0))
    }
  }

  # smooths label 
  smooths<-paste("s(",paste(covs,collapse=","),")", sep="")

  # result
  l<-length(quant)
	Z<- do.call(cbind,lapply(lapply( lapply(1:l, function(x){ 
	  	(diag(ncol(H))-H) %*% getK(datacovs, covs, quant[x], metric, kernel, scale)
	  	} ), function(x){ x * (1/sqrt(sum(diag(t(x) %*% x)))) }) ,unlist))
	colnames(Z)<-1:(l*n)
  lvector<-lapply(1:l,function(x){(((x-1)*n)+1):(x*n)})
  names(lvector)<-paste("quant",quant)
  if (!robust) res<-gt(response=response, alternative=Z, null=null, data=data, ... ,subsets=lvector)
  else res<-gt(response=response, alternative=Z, null=null, data=data, ... )
	res@structure <- list(data.frame(smooths,quant, metric, kernel))
  if (!returnZ) return(res)
  else return(list(result=res, Z=Z))
}

###################################################
# External functions for gtKS
###################################################

getK <- function(data, covs, quant=.25, 
			    metric = c("euclidean", "pearson"), 
				kernel=c("uniform", "exponential", "triangular","neighbours","gauss"), 
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
  if ( kernel == "gauss" ) {
  	dDD <- DD[diag(nrow(DD)) == 0]
    h <- quantile(dDD, quant)
    U <- apply(DD, c(1,2), function(x) { as.numeric( dnorm(x/h) ) } )
  }
  U <- sweep(U, 2, colSums(U), "/")
  U
}

