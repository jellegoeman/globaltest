\name{vandeVijver}

\docType{data}

\alias{vandeVijver}

\title{A subset of the Breast cancer data from Rosetta Inpharmatics LLC and
Netherlands Cancer Institute.}

\description{An \code{\link[Biobase:ExpressionSet-class]{ExpressionSet}}
containing expression data of 100 breast cancer patients on 407
selected genes plus the following phenotype data.
  \describe{
  \item{SampleID}{Sample ID given by NKI for 295 samples in cohort
study (NEJM 347, 01999, 2002)}
  \item{FirstSeriesID}{Sample ID in first series given by NKI for 78 samples
in our first study (Nature 415, p530, 2002)}
  \item{Posnodes}{Lymph node status from pathology report}
  \item{EVENTmeta}{Distant metastasis at any time during follow-up, 0 = no,
1 = yes}
  \item{EVENTdeath}{Death, 0 = no, 1 = yes}
  \item{TIMEsurvival}{Total survival interval in years between first date of
treatment and last date of follow-up (if EVENTdeath = 0) or date of
death (if EVENTdeath = 1)}
  \item{TIMErecurrence}{Disease free interval in years between first date
of treatment and last date of follow-up (if all 'EVENTs' are 0) or
date of first EVENT (if one ore more 'EVENTs' are 1, i.e., local
recurrence, loco-regional recurrence, second primary tumor or
distant metastasis)}
  \item{TIMEmeta}{Metastasis free interval in years between first date of
treatment and date of diagnosis of distant metastasis (only if
EVENTmeta = 1) (note that interval relates to metastasis at any time
during follow-up and not just as first event)}
  \item{ESR1}{Estrogen receptor alpha expression measurement from
microarray (0 = no, 1 = yes)}
  \item{NIH}{Patient characteristics according to NIH criteria, 1 = low
risk, 0 = high risk}
  \item{St Gallen}{Patient characteristics according to St. Gallen
criteria, 1 = low risk, 0 = intermediate/high risk}
  }
}

\usage{data(vandeVijver)}

\seealso{\code{\link{annotation.vandeVijver}}.}

\format{An \code{\link[Biobase:ExpressionSet-class]{ExpressionSet}}. See the Biobase package}

\references{M.J. van de Vijver, Y.D. He, L.J. van 't Veer, H. Dai,
A.A.M. Hart, D.W. Voskuil, G.J. Schreiber, J.L. Peterse, C. Roberts,
M.J. Marton, M. Parrish, D. Atsma, A. Witteveen, A. Glas, L.
Delahaye, T. van der Velde, H. Bartelink, S. Rodenhuis, E.T.
Rutgers, S.H. Friend, and R. Bernards (2002). A gene-expression
signature as a predictor of survival in breast cancer. New England
Journal of Medicine 347 (25), 1999--2009.}

\keyword{datasets}
