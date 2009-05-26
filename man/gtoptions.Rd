\name{gt.options}

\alias{gt.options}

\title{Options for globaltest package}

\description{Sets various global options for the functions in the globaltest package.}

\usage{
gt.options (trace, trim, transpose, max.print, warn.deprecated)
}

\arguments{
    \item{trace}{(Default: \code{interactive()}). If \code{TRUE}, prints progress information whenever many tests are to be performed. Such printing of progress information involves the printing of backspace characters, which is not compatible with use of \code{Sweave}.}
    \item{trim}{(Default: \code{FALSE}). If \code{FALSE}, returns an error if covariates in the \code{subsets} argument of \code{gt} are not present in the data; if \code{TRUE}, silently removes these covariates, and remove duplicates.}
    \item{transpose}{(Default: \code{FALSE}). If \code{TRUE}, \code{gt} expects the transposed data format that is usual in genomics, in which the subjects correspond to the columns of the data matrix and the covariates (or probes) are the rows. If \code{FALSE}, \code{gt} expects the usual statistical format instead.}
    \item{max.print}{(Default: \code{Inf}). The maximum number of characters to print for the alias.}
    \item{warn.deprecated}{(Default: \code{TRUE}). Whether or not to give a warning when the deprecated \code{globaltest} function is called.}
}

\details{The globaltest options can be set during a session, and apply to all calls to functions in the globaltest package for the rest of the session. They are not remembered between sessions.
}

\author{Jelle Goeman: \email{j.j.goeman@lumc.nl}}

\seealso{The \code{\link{gt}} function.}

\examples{
  # setting options
  gt.options(max.print=45, trim=TRUE)

  # reading options
  gt.options()
}

\keyword{htest}
