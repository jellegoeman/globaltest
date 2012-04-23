\name{Diagnostic plots for globaltest}

\alias{covariates}
\alias{features}
\alias{subjects}
\alias{diagnostics}


\title{Global Test diagnostic plots}

\description{Plots to visualize the result of a Global Test in terms of the contributions of the covariates and the subjects.}

\usage{
covariates(object,
            what = c("p-value", "statistic", "z-score", "weighted"),
            cluster = "average", alpha = 0.05, sort = TRUE, zoom = FALSE,
            legend = TRUE, plot = TRUE, colors, alias, help.lines = FALSE,
            cex.labels = 0.6, pdf, trace)

features(...)

subjects(object,
            what = c("p-value", "statistic", "z-score", "weighted"),
            cluster = "average", sort = TRUE, mirror = TRUE,
            legend = TRUE, colors, alias, help.lines = FALSE,
            cex.labels = 0.6, pdf)
}


\arguments{
    \item{object}{A \code{\link{gt.object}}, usually created by a call to \code{\link{gt}}. The object must contain only a single test result, unless the \code{pdf} argument is used. See the help page of \code{\link{gt.object}} on reducing such an object in case it contains more than one test.}
    \item{what}{Gives a choice between various presentations of the same plot. See below under details.}
    \item{cluster}{The type of hierarchical clustering performed for the dendrogram. Default is average linkage clustering. For other options, see \code{\link[stats:hclust]{hclust}}. Setting \code{cluster = "none"} or \code{cluster = FALSE} suppresses the dendrogram altogether.}
    \item{alpha}{Parameter between 0 and 1. Sets the level of family-wise error control in the multiple testing procedure performed on the dendrogram. See below under details.}
    \item{sort}{If \code{TRUE}, the plot sorts the bars with the most significant covariates and subjects to the left, as far as is possible within the constraints of the dendrogram (if present).}
    \item{zoom}{If \code{TRUE}, discards non-significant branches from the dendrogram with the corresponding covariates. This is especially useful for large sets to "zoom" in on the significant results. If no dendrogram is requested, \code{zoom = TRUE} discards all covariates that are not significant after Holm multiple testing correction.}
    \item{legend}{If \code{TRUE}, draws a legend in the plot. To override the default labels of the legend, \code{legend} may also be given as a character vector with the labels of the legend.}
    \item{plot}{If \code{FALSE}, suppress all plotting.}
    \item{colors}{The colors to be used for the bars. See \code{\link[grDevices:rgb]{rgb}} for details on color specification.}
    \item{alias}{Optional alternative labels for the bars in the plots. Should be a character vector of the same length as the number of covariates or subjects, respectively.}
    \item{help.lines}{If \code{TRUE}, prints grey dotted lines that help connect the dendrogram to the bars.}
    \item{cex.labels}{Magnification factor for the x-axis labels.}
    \item{pdf}{Optional filename (\code{character}) of the pdf file to which the plots are to be written. If a filename is provided in \code{pdf}, many \code{covariates} or \code{subjects} plots of multiple tests can be made with a single call to \code{covariates} or \code{subjects}, writing the results to a pdf file.}
    \item{trace}{If \code{TRUE}, prints progress information. Note that printing progress information involves printing of backspace characters, which is not compatible with use of Sweave. Defaults to \code{\link{gt.options}}\code{()$trace}.}
    \item{mirror}{If \code{TRUE}, plots the reverse of the scores for the subjects with negative residual response, so that "good" scores are positive for all subjects.}
    \item{...}{All arguments of \code{features} are identical to those of \code{covariates}.}
}

\details{These two diagnostic plots decompose the test statistics into the contributions of the different covariates and subjects to make the influence of these covariates and subjects visible.

The \code{covariates} plot exploits the fact that the global test statistic for a set of alternative covariates can be written as a weighted sum of the global test statistics for each single contributing covariate. By displaying these component global test results in a bar plot the \code{covariates} plot gives insight into the subset of covariates that is most responsible for the significant test result. The plot can show the \code{p-values} of the component tests on a reversed log scale (the default); their test \code{statistics}, with stripes showing their mean and standard deviation under the null hypothesis; the \code{z-scores} of these test statistics, standardized to mean zero and standard deviation one; or the \code{weighted} test statistics, where the test statistics are multiplied by the relative weight that each covariate carries in the overall test. See the Vignette for more details.

The dendrogram of the \code{covariates} plot is based on correlation distance if the \code{directional} argument was set to \code{TRUE} in the call to \code{\link{gt}}, and uses absolute correlation distance otherwise. The coloring of the dendrogram is based on the multiple testing procedure of Meinshausen (2008): this procedure controls the family-wise error rate on all \code{2n-1} hypotheses associated with the subsets of covariates induced by the clustering graph. All significant subsets are colored black; non-significant ones remain grey. This coloring serves as an additional aid to find the subsets of the covariates most contributing to a significant test result.

The \code{features} function is a synonym for \code{covariates}, using exactly the same arguments.

The \code{subjects} plot exploits the fact that the global test can be written as a sum of contributions of each individual. Each of these contributions is itself a test statistic for the same null hypothesis as the full global test, but one which puts a greater weight on the observed information of a specific subject. These test statistic of subject \code{i} is significant if, for the other subjects, similarity in the alternative covariates to subject \code{i} tends to coincide with similarity in residual response to subject \code{i}. Like the \code{covariates} plot, the \code{subjects} plot can show the \code{p-values} of these component tests on a reversed log scale (the default); their test \code{statistics}, with stripes showing their mean and standard deviation under the null hypothesis; the \code{z-scores} of these test statistics, standardized to mean zero and standard deviation one; or the \code{weighted} test statistics, where the test statistics are multiplied by the relative weight that each covariate carries in the overall test. Setting \code{mirror=FALSE} reverses the bars of subjects with a negative residual response (not applicable if \code{p-values} are plotted). The resulting \code{statistics} values have the additional interpretation that they are proportional to the first order estimates of the linear predictors of each subject under the alternative, i.e. subjects with positive values have higher means under the alternative than under the null, and subjects with negative values have lower means under the alternative than under the null. See the Vignette for more details.

The dendrogram of the \code{subjects} plot is always based on correlation distance. There is no analogue to Meinshausen's multiple testing method for this dendrogram, so multiple testing is not performed.
}

\note{The term "z-score" is not meant to imply a normal distribution, but just refers to a studentized score. The z-scores of the \code{subjects} plot are asymptotically normal under the null hypothesis; the z-scores of the \code{covariates} plot are asymptotically distributed as a chi-squared variable with one degree of freedom.
}

\value{If called to make a single plot, the \code{covariates} function returns an object of class \code{\link{gt.object}}. Several methods are available to access this object: see \code{\link{gt.object}}. The \code{subjects} function returns a matrix. If called to make multiple plots, both functions return \code{NULL}.}

\references{
General theory and properties of the global test are described in

Goeman, Van de Geer and Van Houwelingen (2006) Journal of the Royal Statistical Society, Series B 68 (3) 477-493.

Meinshausen's method for multiple testing

Meinshausen (2008) Biometrika 95 (2) 265-278.

For more references related to applications of the test, see the vignette GlobalTest.pdf included with this package.}

\author{Jelle Goeman: \email{j.j.goeman@lumc.nl}; Livio Finos.}

\seealso{
Diagnostic plots: \code{\link{covariates}}, \code{\link{subjects}}.

The \code{\link{gt.object}} function and useful functions associated with that object.

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

    # Preparation: test
    res <- gt(Y,X)

    # Covariates
    covariates(res)
    covariates(res, what = "w")
    covariates(res, zoom = TRUE)

    # Subjects
    subjects(res)
    subjects(res, what = "w", mirror = FALSE)

    # Change legend, colors or labels
    covariates(res, legend = c("upregulated", "downregulated"))
    covariates(res, col = rainbow(2))
    covariates(res, alias = letters[1:10])

    # Extract data from the plot
    out <- covariates(res)
    result(out)
    extract(out)
}

\keyword{htest}
