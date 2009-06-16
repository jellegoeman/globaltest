\name{gt}

\alias{gt}

\title{Global Test (new implementation, 2009)}

\description{Tests a low-dimensional null hypothesis against a potentially high-dimensional alternative in regression models (linear regression, logistic regression, poisson regression, Cox proportional hazards model).}

\usage{
gt(response, alternative, null, data, test.value, 
     model = c("linear", "logistic", "cox", "poisson", "multinomial"), levels,
     directional = FALSE, standardize = FALSE, permutations = 0, subsets, 
     weights, alias, x = FALSE, trace)
}

\arguments{
    \item{response}{The response vector of the regression model. May be supplied as a vector or as a \code{\link[stats:formula]{formula}} object. In the latter case, the right hand side of \code{response} is passed on to \code{alternative} if that argument is missing, or otherwise to \code{null}.}
    \item{alternative}{The part of the design matrix corresponding to the alternative hypothesis. The covariates of the null model do not have to be supplied again here. May be given as a half \code{\link[stats:formula]{formula}} object (e.g. \code{~a+b}). In that case the intercept is always suppressed. Alternatively, the \code{alternative} argument may also be given as an \code{\link[Biobase:ExpressionSet]{ExpressionSet}} object, in which case \code{t(exprs(alternative))} is used for the \code{alternative} argument, and \code{pData(alternative)} is passed on to the \code{data} argument if that argument is missing.}
    \item{null}{The part of the design matrix corresponding to the null hypothesis. May be given as a design matrix or as a half \code{\link[stats:formula]{formula}} object (e.g. \code{~a+b}). The default for \code{null} is \code{~1}, i.e. only an intercept. This intercept may be suppressed, if desired, with \code{null = ~0}.}
    \item{data}{Only used when \code{response}, \code{alternative}, or \code{null} is given in formula form. An optional data frame, list or environment containing the variables used in the formulae. If the variables in a formula are not found in \code{data}, the variables are taken from environment(formula), typically the environment from which \code{gt} is called.}
    \item{test.value}{An optional vector regression coefficients to test. The default is to test the null hypothesis that all regression coefficients of the covariates of the alternative are zero. The \code{test.value} argument can be used to test a value other than zero. The coefficients are applied to the design matrix of \code{alternative} before any standardization (see the \code{standardize} argument).}
    \item{model}{The type of regression model to be tested. If omitted, the function will try to determine the model from the class and values of the \code{response} argument.}
    \item{levels}{Only used if response is \code{\link{factor}}. Selects a subset of \code{levels(response)} to be tested, given as a character vector. If a vector of length \code{>1}, the test uses only the subjects with the specified outcome categories. If \code{levels} is of length 1, the test reduces the response to a two-valued factor, testing the specified outcome category against the others combined.}
    \item{directional}{If set to \code{TRUE}, directs the power of the test especially against the alternative that the true regression coefficients under the alternative have the same sign. The default is that the power of the test does not depend on the sign of the true regression coefficients. Set negative \code{weights} for covariates that are expected to have opposite sign.}
    \item{standardize}{If set to \code{TRUE}, standardizes all covariates of the alternative to have unit second central moment. This makes sure that the test result is independent of the relative scaling of the covariates. The default is to let covariates with more variance have a greater weight in the test.}
    \item{permutations}{The number of permutations to use. The default, \code{permutations = 0}, uses the asymptotic distribution. The asymptotic distribution is the exact distribution in case of the linear model with normal errors.}
    \item{subsets}{Optional argument that can be used to test one or more subsets of the covariates in \code{alternative}. Can be a vector of column names or column indices of \code{alternative}, or a list of such vectors. In the latter case, a separate test will be performed for each subset.}
    \item{weights}{Optional argument that can be used to give certain covariates in \code{alternative} greater weight in the test. Can be a vector or a list of vectors. In the latter case, a separate test will be performed for each weight vector. If both \code{subsets} and \code{weights} are specified as a list, they must have the same length. In that case, \code{weights} vectors may have either the same length as the number of covariates in \code{alternative}, or the same length as the corresponding subset vector. Weights can be negative; the sign has no effect unless \code{directional} is \code{TRUE}.}
    \item{alias}{Optional second label for each test. Should be a vector of the same length as \code{subsets}. See also \code{\link{alias}}.}
    \item{x}{If \code{TRUE}, gives back the \code{null} and \code{alternative} design matrices. Default is not to return these matrices.}
    \item{trace}{If \code{TRUE}, prints progress information. This is useful if many tests are performed, i.e.\ if \code{subsets} or \code{weights} is a list. Note that printing progress information involves printing of backspace characters, which is not compatible with use of Sweave. Defaults to \code{\link{gt.options}}\code{()$trace}.}
}

\details{The Global Test tests a low-dimensional null hypothesis against a (potentially) high-dimensional alternative, using the locally most powerful test of Goeman et al (2006). In this regression model implementation, it tests the null hypothesis \code{response ~ null}, that the covariates in \code{alternative} are not associated with the response, against the alternative model \code{response ~ null + alternative} that they are.

The test has a wide range of applications. In \link[gtGO]{gene set testing} in microarray data analysis \code{alternative} may be a matrix of gene expression measurements, and the aim is to find which of a collection of predefined \code{subsets} of the genes (e.g. Gene Ontology terms or KEGG pathways) is most associated with the \code{response}. In penalized regression or other machine learning techniques, \code{alternative} may be a collection of predictor variables that may be used to predict a \code{response}, and the test may function as a useful pre-test to see if training the classifier is worthwhile. In goodness-of-fit testing, \code{null} may be a model with linear terms fitted to the \code{response}, and \code{alternative} may be a large collection of non-linear terms. The test may be used in this case to test the fit of the null model with linear terms against a non-linear alternative.

See the vignette for extensive examples of these applications.
}

\note{If \code{null} is supplied as a \code{\link[stats:formula]{formula}} object, an intercept is automatically included. As a consequence \code{gt(Y, X, Z)} will usually give a different result from \code{gt(Y, X, ~Z)}. The first call is equivalent to \code{gt(Y, X, ~0+Z)}, whereas the second call is equivalent to \code{gt(Y, X, cbind(1,Z))}.

P-values from the asymptotic distribution are accurate to at least two decimal places up to a value of around \code{1e-12}. Lower p-values are numerically less reliable.

Missing values are allowed in the \code{alternative} matrix only. Missing values are imputed conservatively (i.e. under the null hypothesis). Covariates with many missing values get reduced variance and therefore automatically carry less weight in the test result.
}

\value{The function returns an object of class \code{\link[gt.object-class]{gt.object}}. Several operations and diagnostic plots can be made from this object. See also \link[diagnostics]{Diagnostic plots}.}

\references{
General theory and properties of the global test are described in

Goeman, Van de Geer and Van Houwelingen (2006) Journal of the Royal Statistical Society, Series B 68 (3) 477-493.

For references related to applications of the test, see the vignette GlobalTest.pdf included with this package.}

\author{Jelle Goeman: \email{j.j.goeman@lumc.nl}; Jan Oosting}

\seealso{
Diagnostic plots: \code{\link{covariates}}, \code{\link{subjects}}.

The \code{\link[gt.object-class]{gt.object}} function and useful functions associated with that object.

Many more examples in the vignette!
}

\examples{
    # Simple examples with random data here
    # Real data examples in the Vignette

    # Random data: covariates A,B,C are correlated with Y
    set.seed(1)
    Y <- rnorm(20)
    X <- matrix(rnorm(200), 20, 10)
    X[,1:3] <- X[,1:3] + Y
    colnames(X) <- LETTERS[1:10]

    # Compare the global test with the F-test
    gt(Y, X)
    anova(lm(Y~X))

    # Using formula input
    res <- gt(Y, ~A+B, null=~C+E, data=data.frame(X))
    summary(res)

    # Beware: null models with and without intercept
    Z <- rnorm(20)
    summary(gt(Y, X, null=~Z))
    summary(gt(Y, X, null=Z))

    # Logistic regression
    gt(Y>0, X)

    # Subsets and weights (1)
    my.sets <- list(c("A", "B"), c("C","D"), c("D", "E"))
    gt(Y, X, subsets = my.sets)
    my.weights <- list(1:2, 2:1, 3:2)
    gt(Y, X, subsets = my.sets, weights=my.weights)

    # Subsets and weights (2)
    gt(Y, X, subset = c("A", "B"))
    gt(Y, X, subset = c("A", "A", "B"))
    gt(Y, X, subset = c("A", "A", "B"), weight = c(.5,.5,1))

    # Permutation testing
    summary(gt(Y, X, perm=1e4))
}

\keyword{htest}
