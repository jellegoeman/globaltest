\name{gt goodness of fit methods}

\alias{gtPS}
\alias{gtK}
\alias{gtI}
\alias{bbase}
\alias{repdes}
\alias{getPS}
\alias{getK}

\title{Goodness of fit testing in regression models using Global Test}

\description{
Tests the goodness of fit of a regression model against a specified alternative using the Global Test. 
Three main functions are provided: \code{gtPS} uses Penalized Splines, \code{gtK} uses Kernel smoothers and \code{gtI} uses linear Interactions. The other functions are for external use in combination with \code{gt}.}

\usage{
gtPS(response, null, data, model, ..., covs, bdeg = 3, nint= 10, pord = 2,
	interact = FALSE, by = NULL, robust = FALSE, termlabels = FALSE)
	
gtK(response, null, data, model, ..., covs, quant = .25, 
       metric = c("euclidean", "pearson"), 
       kernel=c("uniform", "exponential", "triangular", "neighbours", "Gauss"),
       robust = FALSE, scale = TRUE, termlabels = FALSE)
	
gtI(response, null, data, model, ..., covs, iorder=2)

bbase(x, bdeg, nint)

repdes(B, pord, K=NULL, tol = 1e-10)

getPS(data, covs, bdeg=3, nint=10, pord=2, by=NULL, interact=FALSE)

getK(data, covs, quant=.25, metric = c("euclidean", "pearson"), 
	kernel=c("uniform", "exponential", "triangular","neighbours","Gauss"), 
	scale=TRUE )

}

\arguments{
    \item{response}{The response vector of the regression model. May be
      supplied as a vector, as a \code{\link[stats:formula]{formula}} object, or as an object of class \code{\link[stats:lm]{lm}}, \code{\link[stats:glm]{glm}} or \code{\link[stats:lm]{coxph}}. In the last two cases, the specification of \code{null} is not required.}
    \item{null}{The null design matrix. May be given as a matrix or as a half \code{\link[stats:formula]{formula}} object (e.g. \code{~a+b}).}
    \item{data}{Only used when \code{response} or \code{null} is given in formula form. An optional data frame, list or environment containing the variables used in the formulae.}
    \item{model}{The type of regression model to be tested. If omitted, the function will try to determine the model from the class and values of the \code{response} argument.}
    \item{...}{Any other arguments are also passed on to \code{\link{gt}}.}
    \item{covs}{A variable or a vector of variables that are the covariates the smooth terms are function of.}
    \item{bdeg}{A vector or a list of vectors which specifies the degree of the B-spline basis.}
    \item{nint}{A vector or a list of vectors which specifies  the number of intervals determined by equally-spaced knots.}
    \item{pord}{A vector or a list of vectors which specifies  the order of the differences indicating the type of the penalty imposed to the coefficients.}
    \item{interact}{\code{TRUE} to consider a multidimensional smooth function of \code{covs}.}
    \item{by}{A factor variable or a vector of factor variables to let the smooth terms specified in \code{covs} interact with factors terms. }
    \item{termlabels}{\code{TRUE} to consider e.g. \code{s(log(cov))} instead of \code{s(cov)} when \code{null=~ log(cov)} and \code{covs} is missing.}
    \item{robust}{\code{TRUE} to obtain an overall test which combines multiple specifications of the B-spline basis arguments (when \code{bdeg}, \code{nint} and \code{pord} are lists) or multiple specifications of the bandwidth (when \code{quant} is a vector of quantiles).}
    \item{quant}{The smoothing bandwidth to be used, expressed as the percentile of the distribution of distance between observations, with default the 25th percentile. To investigate the sensitivity to different choices, \code{quant} can be a vector of percentiles. See also \code{robust} argument.}
    \item{metric}{A character string specifying the metric to be used. The available options are "euclidean" (the default), "pearson" and "mixed" (to be implemented).
"mixed" distance is chosen automatically if some of the selected covariates are not numeric.}
    \item{kernel}{A character string giving the smoothing kernel to be used. This must be one of "uniform", "exponential", "triangular", "neighbours", or "Gauss", with default "uniform".}
    \item{scale}{\code{TRUE} to scale the covariates before computing the distance.}
    \item{iorder}{Order of the linear interactions, e.g. second order interactions, third order etc.}
    \item{x}{A numeric vector of values at which to evaluate the B-spline basis}
    \item{B}{Matrix containing the B-spline basis}
    \item{K}{Penalty matrix (i.e. the penalty term is the quadratic form of K and the spline coefficients)}
    \item{tol}{Eigenvalues smaller than \code{tol} are considered zero}
}
\details{


These are functions to test for specific types of lack of fit by using the Global Test. 
Suppose that we are concerned with the adequacy of some regression model \code{response ~ null}, such as \code{Y ~ X1 + X2}.
The alternative model can be cast into the generic form \code{response ~ null + alternative}, which comprises different models corresponding to different types of lack of fit. Thus, the specification of \code{alternative} is required. It identifies the type of lack of fit the test is directed against.

By using \code{gtPS}, the alternative is given by a user specified sum of smooth functions of continuous covariates, 
e.g.  \code{alternative= ~ s(X1)} when \code{covs="X1"} and \code{alternative= ~ s(X1) + s(X2)} when \code{covs=c("X1","X2")}.
Smooth terms are represented using P-splines as proposed by Eilers and Marx (1996). This approach consists in constructing 
a B-spline basis of degree \code{bdeg} with \code{nint + 1} equidistant knots, where a difference penalty of order \code{pord} is applied to the basis coefficients. Any sane combination of penalty order and basis degree is allowed. 
If \code{interact=TRUE}, the alternative is given by a multidimensional smooth function of \code{covs}, which is represented by a tensor product of the marginal B-splines bases and the Kronecker sum of the marginal penalties, e.g. \code{alternative= ~ s(X1,X2)} when \code{covs=c("X1","X2")} and \code{interact=TRUE}. If the \code{null} model includes continuous and factor variables, e.g. \code{X1} continuous and \code{X2} factor, 
the argument \code{by} allows to construct varying-coefficient alternative models, e.g. 
\code{alternative= ~ s(X1):X2}, when \code{covs="X1"} and \code{by="X2"}. If the response is survival and the \code{null} model includes factor variables, e.g. \code{response=Surv(time,status)} and \code{X2} is a factor, the argument \code{by} allows to construct time-varying effects alternative models, e.g. \code{alternative= ~ s(time):X2} by using \code{covs="time"} and \code{by="X2"}. 

By using \code{gtK} the alternative is given by a user specified multidimensional smooth term, e.g. \code{alternative= ~ s(X1, X2)} when \code{covs=c("X1","X2")}. 
Multidimensional smooth terms are represented by a kernel smoother defined by a distance measure (\code{metric}), a kernel shape (\code{kernel}) and a bandwidth (\code{quant}). 
Because the test is sensitive to the chosen value of \code{quant}, it is possible to specify \code{quant} as a vector of different values in combination with \code{robust=TRUE}. 
Distance measures for factor covariates and for the situation that both continuous and factor covariates are present are constructed as in le Cessie and van Houwelingen (1995), e.g. \code{covs=c("X1","X2")} and \code{distance="mixed"} when \code{X1} continuous and \code{X2} factor. 

By using \code{gtI}, the alternative is given by all the possible ith-order linear interactions between \code{covs}, e.g. \code{alternative= ~ X1:X2} when \code{covs=c("X1","X2")} and \code{iorder=2}. 

The remaining functions are meant for constructing the alternative design matrix that will be used in the \code{alternative} argument of the \code{interact=TRUE} function. 
\code{bbase} constructs the B-spline basis for the covariate \code{cov} and \code{ndiff} builts a difference penalty matrix of order \code{pord}. 
These function are based on the functions provided by Eilers and Marx (1996). 
\code{repdes} reparameterizes a spline basis \code{B} via the spectral decomposition of its penalty matrix \code{K} into a design for
the unpenalized part of the function and a design for the penalized part with i.i.d. coefficients. 
\code{getPS} and \code{getK} return a reparameterized B-spline basis and a Kernel smoothing matrix, respectively. 

See the vignette for more examples.
}


\note{Currently linear (normal), logistic, multinomial logistic and Poisson regression models with canonical links and Cox's proportional hazards regression model are supported.
}

\value{The function returns an object of class \code{\link{gt.object}}. Several operations and diagnostic plots can be made from this object.}

\references{

Eilers, Marx (1996).
Flexible smoothing with B-splines and penalties.
Statistical Science, 11: 89-121.

le Cessie, van Houwelingen (1995). Testing the Fit of a Regression Model Via Score Tests in Random Effects Models.
Biometrics 51: 600-614.


For references related to applications of the test, see the vignette GlobalTest.pdf included with this package.}

\author{Aldo Solari: \email{aldo.solari@unimib.it}}

\seealso{
The \code{\link{gt}} function. The \code{\link{gt.object}} and useful functions associated with that object.
}

\examples{

    # Random data: univariate 
    set.seed(0)
    X1<-runif(50)
    s1 <- function(x) exp(2 * x)
    e <- rnorm(50)
    Y <-  s1(X1) + e
    
    ### gtPS
    gtPS(Y~X1)
    
    # model input
    rdata<-data.frame(Y,X1)
    nullmod<-lm(Y~X1,data=rdata)
    gtPS(nullmod)
    
    # formula input and termlabels
    gtPS(Y~exp(X1),data=rdata)
    gtPS(Y~exp(X1),covs="exp(X1)",data=rdata)
    gtPS(Y~exp(X1),data=rdata, termlabels=TRUE)
    
    # P-splines arguments
    gtPS(Y~X1, nint=list(a=2, b=10, c=50))
    gtPS(Y~X1, nint=list(a=2, b=10, c=50), pord=0)
    gtPS(Y~X1, nint=list(a=10, b=10), pord=list(a=0,b=2))
    gtPS(Y~X1, nint=list(a=10, b=10), pord=list(a=0,b=2), robust=TRUE)
        
    # Random data: additive model 
    X2<-runif(50)
    s2 <- function(x) 0.2 * x^11 * (10 * (1 - x))^6 + 10 * (10 * x)^3 * (1 - x)^10
    Y <-  s1(X1) + s2(X2) + e
    gtPS(Y~X1+X2)
    gtPS(Y~X1+X2, pord=list(a=c(0,2), b=c(2,0)))
    gtPS(Y~X1+X2, covs="X2")
    
    ### gtK  
    gtK(Y~X1+X2)
    gtK(Y~X1+X2, quant=seq(.05,.95,.1))
    gtK(Y~X1+X2, quant=seq(.05,.95,.1), robust=TRUE)
    
    # Random data: smooth surface
    s12 <- function(a, b, sa = 1, sb = 1) {
            (pi^sa * sb) * (1.2 * exp(-(a - 0.2)^2/sa^2 - (b - 0.3)^2/sb^2) + 
            0.8 * exp(-(a - 0.7)^2/sa^2 - (b - 0.8)^2/sb^2))
            }
    Y <- s12(X1,X2) + e
    
    # Non-linear interaction
    gtPS(Y~X1*X2)
    gtPS(Y~X1*X2, interact=TRUE)
    gtK(Y~X1*X2, quant=seq(.05,.95,.1), robust=TRUE)
        
    # Time-varying effects
    require("survival")
    data(rats)
    gtPS(Surv(time,status)~rx,covs="time", by="rx",data=rats, model="cox", pord=1)

}

\keyword{goodness of fit}











