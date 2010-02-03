\name{Multiple testing methods of the globaltest package}

\alias{focusLevel}
\alias{findFocus}
\alias{inheritance}
\alias{leafNodes}
\alias{draw}


\title{Multiple testing correction for the Global Test}

\description{A collection of multiple testing procedures for the Global Test. Methods for the focus level procedure of Goeman and Mansmann for graph-structured hypotheses, and for the inheritance procedure based on Meinshausen.}

\usage{
# The focus level method:
focusLevel (test, sets, focus, ancestors, offspring,
           stop = 1, atoms = TRUE, trace)

findFocus (sets, ancestors, offspring, maxsize = 10, atoms = TRUE)


# The inheritance method:
inheritance (test, sets, weights, ancestors, offspring, Shaffer,
            homogeneous = FALSE, trace)


# Utilities for focus level and inheritance method:
leafNodes (object, alpha=0.05, type = c("focuslevel","inheritance"))

draw (object, alpha = 0.05, type = c("focuslevel","inheritance"),
      names=FALSE, sign.only = FALSE, interactive = FALSE)
}

\arguments{
    \item{object}{A \code{\link{gt.object}}, usually one in which more than one test was performed.}
    \item{test}{Either a \code{\link{function}} or \code{\link{gt.object}}. If a function, that function should take as its argument a vector of covariate labels, and return  (raw) p-value. See the examples below. If a \code{\link{gt.object}} the call to \code{\link{gt}} that created it must have had all the covariates of \code{sets} (below) in its \code{alternative} argument.}
    \item{sets}{A named \code{\link[base:list]{list}} representing covariate sets of the hypotheses of interest, for which adjusted p-values are to be calculated. If it is missing but \code{test} is
        a \code{\link{gt.object}}, the \code{subsets} slot of that
	object will be used. If used in the \code{inheritance},
	\code{sets} describe a tree structure of hypotheses. In this
	case, object of class \code{\link[stats:hclust]{hclust}} or
	\code{\link[stats:dendrogram]{dendrogram}}.}

    \item{focus}{The focus level of the focus level method. Must be a subset of \code{names(sets)}. Represents the level of the graph at which the method is focused,
        i.e. has most power.}
    \item{ancestors}{An environment or list that maps each set in \code{sets} to all its ancestors, i.e. its proper supersets. If missing, \code{ancestors} is determined
        from the input of \code{offspring}, or, if that is also missing, from the input of \code{sets} (time-consuming).}
    \item{offspring}{An environment or list that maps each set in \code{sets} to all its offspring, i.e. its proper subsets. If missing, \code{offspring} is determined
        from the input of \code{ancestors}, or, if that is also missing, from the input of \code{sets} (time-consuming).}
    \item{stop}{Determines when to stop the algorithm. If \code{stop} is set to a value smaller than or equal to 1, the algorithm only calculates familywise error rate corrected
        p-values of at most \code{stop}. If \code{stop} is set to a value greater than 1, the algorithm stops when it has rejected at least \code{stop} hypotheses. If set
        to exactly 1, the algorithm calculates all familywise error rate corrected p-values. Corrected p-values that are not calculated are reported as \code{NA}.}
    \item{atoms}{If set to \code{TRUE}, the focus level algorithm partitions the offspring of each focus level set into the smallest possible building blocks, called
        atoms. Doing this often greatly accelerates computation, but sometimes at the cost of some power.}
    \item{trace}{If set to \code{TRUE}, reports progress information. The default is obtained from \code{\link[=gt.options]{gt.options()$trace}}. Alternatively, 
        setting \code{trace = 2} gives much more extensive output (\code{focusLevel} only).}
    \item{maxsize}{Parameter to choose the height of the focus level. The focus level sets are chosen in such a way that the number of tests that is to be done for each
        focus level set is at most \code{2^maxsize - 1}.}
    \item{alpha}{The alpha level of familywise error control for the significant subgraph.}
    \item{Shaffer}{If set to \code{TRUE}, it applys the Shaffer improvement. If \code{Shaffer} is \code{NULL} and \code{object} is a \code{\link{gt.object}} the procedure checks whether  \code{Shaffer=TRUE} is valid, and sets the value accordingly.}
    \item{weights}{Optional weights vector for the leaves nodes. If it is missing but \code{test} is a \code{\link{gt.object}}, the result of \code{\link{weights}(object)} will be used. In all other cases \code{weights} is set to be uniform among all leaf nodes.}
    \item{homogeneous}{If set to \code{TRUE}, redistributes the alpha of rejected leaf node hypotheses homogeneously over the hypotheses under test, rather than to closest related hypotheses.}
    \item{type}{Argument for specifying which multiple testing correction method should be used. Only relevant if both the inheritance and the focuslevel procedures were
        performed on the same set of test results.}
    \item{names}{If set to \code{TRUE}, draws the graph with node names rather than numbers.}
    \item{sign.only}{If set to \code{TRUE}, draws only the subgraph corresponding to the significant nodes. If \code{FALSE}, draws the full graph with the non-significant
        nodes grayed out.}
    \item{interactive}{If set to \code{TRUE}, creates an interactive graph in which the user can see the node label by clicking on the node.}
}

\details{Multiple testing correction becomes important if the Global Test is performed on many covariate subsets.

If the hypotheses are structured in such a way that many of the tested subsets are subsets of other sets, more powerful procedures can be applied that take advantage of this structure to gain power. Two methods are implemented in the \code{globaltest} package: the \code{inheritance} method for tree-structured hypotheses and the \code{focusLevel} method for general directed acyclic graphs. For simple multiple testing that does not use such structure, see \code{\link[=gt.object]{p.adjust}}.

The \code{focusLevel} procedure makes use of the fact that some sets are subsets or supersets of each other, as specified by the user in the \code{offspring} and
\code{ancestors} arguments. Viewing the subset and superset structure as a graph, the procedure starts testing at a \code{focus} level: a subset of the nodes of the graph.
If the procedure finds significance at this focus level, it proceeds to find significant subsets and supersets of the focus level sets. Like Holm's procedure, the focus
level procedure is valid regardless of the correlation structure between the test statistics.

The focus level method requires the choice of a ``focus level'' in the graph. The \code{findFocus} function is a utility function for automatically choosing a focus level. It chooses a collection of focus level sets in such a way that the number of tests to be done for each focus level node is at most \code{2^maxsize}. In practice this usually means that each focus level node has at most \code{maxsize} leaf nodes as offspring. Choosing focus level nodes with too many offspring nodes may result in excessively long computation times.

The \code{inheritance} method is an alternative method for calculating familywise error rate corrected p-values. Like the focus level method, \code{inheritance} also makes use of the structure of the tested sets to gain power. In this case, however, the graph is restricted to a tree, as can be obtained for example if the tested subsets are obtained from a hierarchical clustering. The inheritance procedure is used in the \code{\link{covariates}} function. Like Holm's method and the focus level method, the inheritance procedure makes no assumptions on the joint distribution of the test statistics.

The \code{leafNodes} function extracts the leaf nodes of the significant subgraph after a focus level procedure was performed. As this graph is defined by its leaf nodes, this is the most efficient summary of the test result. Only implemented for \code{\link{gt.object}} input.

The \code{draw} function draws the graph, displaying the significant nodes. It either draws the full graph with the non-significant nodes grayed out (\code{sign.only = TRUE}), or it draws only the subgraph corresponding to the significant nodes.

See the vignette for extensive applications. }

\note{
In the graph terminology of the focus level method, \code{ancestor} means superset, and \code{offspring} means subset.

The validity of the focus level procedure depends on certain assumptions on the null hypothesis that is tested for each set. See the paper by Goeman and Mansmann (cited below) for the precise assumptions. Similar assumptions are necessary for the Shaffer improvement of the inheritance procedure.}

\value{The function \code{multtest} returns an object of class
  \code{\link[globaltest:gt.object-class]{gt.object}} with an
  appropriate column added to the test results matrix.


The \code{focusLevel} and \code{inheritance} functions returns a
\code{\link{gt.object}} if a
\code{\link[globaltest:gt.object-class]{gt.object}} argument was given
as input, otherwise it returns a matrix with a column of raw p-values
and a column of corrected p-values.


The function \code{leafNodes} returns a
\code{\link[globaltest:gt.object-class]{gt.object}}.


\code{findFocus} returns a character vector.}

\references{ The methods used by multtest:

Holm (1979) A simple sequentially rejective multiple test procedure. Scandinavian Journal of Statistics 6: 65-70.

Benjamini and Hochberg (1995) Controlling the false discovery rate: a practical and powerful approach to multiple testing. Journal of the Royal Statistical Society, Series
B 57: 289-300.

Benjamini and Yekutieli (2001) The control of the false discovery rate in multiple testing under dependency. Annals of Statistics 29 (4) 1165-1188.

The focus level method:

Goeman and Mansmann (2008) Multiple testing on the directed acyclic graph of gene ontology. Bioinformatics 24 (4) 537-544.

The inheritance method:

Meinshausen (2008) Hierarchical testing of variable importance. Biometrika 95 (2), 265-278.

For references related to applications of the test, see the vignette GlobalTest.pdf included with this package.}


\author{Jelle Goeman: \email{j.j.goeman@lumc.nl}; Livio Finos}

\seealso{ The \code{\link{gt}} function. The \code{\link{gt.object}} function and useful functions associated with that object.

Many more examples in the vignette! }

\examples{
    # Simple examples with random data here
    # Real data examples in the Vignette

    # Random data: covariates A,B,C are correlated with Y
    set.seed(1)
    Y <- rnorm(20)
    X <- matrix(rnorm(200), 20, 10)
    X[,1:3] <- X[,1:3] + Y
    colnames(X) <- LETTERS[1:10]

    # Some subsets of interest
    my.sets1 <- list(abc = LETTERS[1:3], cde  = LETTERS[3:5],
                     fgh = LETTERS[6:8], hij = LETTERS[8:10])
    res <- gt(Y, X, subsets = my.sets1)

    # Simple multiple testing
    p.adjust(res)
    p.adjust(res, "BH")

    # A whole structure of sets
    my.sets2 <- as.list(LETTERS[1:10])
    names(my.sets2) <- letters[1:10]
    my.sets3 <- list(all = LETTERS[1:10])
    my.sets <- c(my.sets2,my.sets1,my.sets3)

    # Do the focus level procedure
    # Choose a focus level by hand
    my.focus <- c("abc","cde","fgh","hij")
    # Or automated
    my.focus <- findFocus(my.sets, maxsize = 8)
    resF <- focusLevel(res, sets = my.sets, focus = my.focus)
    leafNodes(resF, alpha = .1)

    # Compare
    p.adjust(resF, "holm")

    # Focus level with a different test
    Ftest <- function(set) anova(lm(Y~X[,set]))[["Pr(>F)"]][1]
    focusLevel(Ftest, sets=my.sets, focus=my.focus)

    # analyze data using inheritance procedure
    res <- gt(Y, X, subsets = list(colnames(X)))
    # define clusters on the covariates X
    hcl=hclust(dist(t(X)))
    # Do inheritance procedure
    resI=inheritance(res, sets = hcl)
    resI
    leafNodes(resI, alpha = .1)

}

\keyword{htest}
