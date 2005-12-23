\name{gt.barplot-class}
\docType{class}
\alias{gt.barplot-class}
\alias{z.score}
\alias{gt.barplot}
\alias{show,gt.barplot-method}
\alias{[,gt.barplot-method}
\alias{length,gt.barplot-method}
\alias{result,gt.barplot-method}
\alias{plot,gt.barplot-method}
\alias{scale,gt.barplot-method}
\alias{z.score,gt.barplot-method}
\alias{names,gt.barplot-method}
\alias{names<-,gt.barplot-method}

\title{Class "gt.barplot" for results of the functions geneplot and sampleplot}

\description{The class gt.result stores the output of a call to
\code{\link{geneplot}} or \code{\link{sampleplot}}.}

\section{Slots}{
  \describe{
    \item{\code{res}:}{Object of class "matrix". The plotted values, with
    a row for each gene or sample, and the following columns: \code{influence}: the height of the bar,
    \code{Einf}: expected influence under the null hypothesis, \code{sd.inf}:
    standard deviation of the influence under the null hypothesis,
    \code{up}: the colour of the plotted bar (1=green; 0=red)}
    \item{\code{drawlabels}:}{Object of class "logical".
    Whether gene or sample labels should be printed at the x-axis.}
    \item{\code{labelsize}:}{Object of class "numeric". The size of the gene/sample labels.}
    \item{\code{legend}:}{The meaning of the red and green colours of the bars.}
    }
}

\section{Methods}{
  \describe{
    \item{show}{(gt.barplot): Summarizes plot in numbers.}
    \item{plot}{(gt.barplot): Redraws the plot. Possible extra arguments:
    \code{labelsize} and \code{drawlabels}.}
    \item{"["}{(gt.barplot): Extracts a subset of the genes or samples.}
    \item{result}{(gt.barplot): Extracts the results matrix as a \code{data.frame}.}
    \item{length}{(gt.barplot): The number of bars (genes or samples) in the object.}
    \item{scale}{(gt.barplot): Scales the bars to unit standard deviation.}
    \item{z.score}{(gt.barplot): Calculates the z.scores for each bar: (influence - expected inf) / sd(inf).}
    \item{names}{(gt.barplot): The names of the genes or samples.}
    \item{names<-}{(gt.result): Changes the names of the genes or samples.}
  }
}

\author{Jelle Goeman: \email{j.j.goeman@lumc.nl}; Jan Oosting}

\seealso{\code{\link{sampleplot}},
\code{\link{geneplot}}, \code{\link{globaltest}}.}

\keyword{methods}
