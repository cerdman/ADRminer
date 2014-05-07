\encoding{UTF-8}
\name{pvPen2}
\alias{pvPen2}
\alias{pvPen2.pvInd}
\title{Penalized regression detection strategy}
\usage{
\method{pvPen2}{pvInd}(object, aeId = "all", covId = NULL, posConst = TRUE, nDrugMax = 20, parallel = require(parallel), nCores = NULL, \dots)
}
\arguments{
  \item{object}{an object of class pvInd}

  \item{aeId}{The label of the adverse event to be
  regressed. By default, the \code{pvPen} will regress all
  ae, one at a time and this can be very time consuming and
  RAM demanding. This parameter can be either the name of
  the ae(s), or the index of the ae column}

  \item{covId}{a character vector indicating which
  covariate have to be used in the analyses.}

  \item{posConst}{If TRUE (default), the regression
  coefficients are constrained to be non non negative, this
  to ensure that the final model only contains drugs
  "increasing" the risk of a given ae.}

  \item{nDrugMax}{Maximum number of drugs to be included in
  the model. In addition to be computationnaly intensive,
  odels with to many drugs are likely to be very unstable.}

  \item{parallel}{Whether parallel computations should be
  used to speed up the calculation. Be careful as it will
  be more RAM demanding. Only available for linux or Mac
  Os}

  \item{nCores}{If parallel is true, allows to specify the
  number of cores to be used for parallel computing.}

  \item{...}{Further arguments to be passed (None for the
  moment)}
}
\description{
pvPen implements a detection strategy based on the use of
the Lasso regression along with BIC type detection
criteria. The methodology follows two main steps: \itemize{
\item Lasso regression (with a positive constraint on the
regression coefficients) for a grid of constraint
parameters \item selection on the drug associated with the
ae of interest according to the BIC.  } This function
heavily relies on the highly efficient \code{glmnet}
package. \code{pvPen} is computationnaly highly demanding.
It is also strongly advised to use a linux or a mac
computer in order to use several cores. Note also that this
regression approach should be used with ae and drugs both
having a reasonable number of reports.
}

