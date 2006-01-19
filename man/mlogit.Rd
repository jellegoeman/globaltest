\name{mlogit}

\alias{mlogit}
\alias{mlogit-class}
\alias{residuals,mlogit-method}
\alias{fitted.values,mlogit-method}
\alias{coefficients,mlogit-method}
\alias{summary,mlogit-method}
\alias{show,mlogit-method}

\title{Multinomial Logistic Regression}

\description{Fits a multinomial logistic regression model to a nominal scale outcome.}

\usage{mlogit(formula, data, control = glm.control())}

\arguments{
    \item{formula}{An object of class \code{\link[stats:formula]{formula}}
containing a symbolic description of the model to be fit. See the documentation
of \code{\link[stats:formula]{formula}} for details.}
    \item{data}{An optional data frame containing the variables in the model.
If not found in 'data', the variables are taken from the environment from which
'mlogit' is called.}
    \item{control}{A list of parameters for controlling the fitting process.
See the documentation of \code{\link[stats:glm.control]{glm.control}} for details.} }

\details{The function mlogit fits a multinomial logistic regression
model for a multi-valued outcome with nominal scale. The
implementation and behaviour are designed to mimic those of
\code{\link[stats:glm]{glm}}, but the options are (as yet) more
limited. Missing values are not allowed in the data.

The model is fitted without using a reference outcome category; the
parameters are made identifiable by the requirement that the sum of
corresponding regression coefficients over the outcome categories is
zero.}

\value{An object of (S4) class \code{mlogit}. The class has slots:
coefficients (matrix), standard.err (matrix), fitted.values
(matrix), x (matrix), y (matrix), formula (formula), call (call),
df.null (numeric), df.residual (numeric), null.deviance (numeric),
deviance (numeric), iter (numeric), converged (logical).

Methods implemented for the \code{mlogit} class are
\code{coefficients}, \code{fitted.values}, \code{residuals} and
which extract the relevant quantities, and \code{summary}, which
gives the same output as with a \code{\link[stats:glm]{glm}}
object.}

\author{Jelle Goeman: \email{j.j.goeman@lumc.nl}; Jan Oosting}

\seealso{\code{\link{glm}}, \code{\link[nnet:multinom]{multinom}}.}

\examples{
  y <- factor(rep(1:4,5))
  x <- 1:20
  fit <- mlogit(y ~ x)
  summary(fit)
  residuals(fit)
}

\keyword{nonlinear}
