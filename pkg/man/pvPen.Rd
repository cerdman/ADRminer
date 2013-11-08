\encoding{UTF-8}
\name{pvPen}
\alias{pvPen}
\alias{pvPen.pvInd}
\title{Penalized regression detection strategy}
\usage{
  \method{pvPen}{pvInd}(object, aeLab="all", cov=NULL,
    detectCriter=c("BIC", "AIC"), posConst=TRUE,
    nDrugMax=20, parallel=FALSE, nCore=NULL, \dots)
}
\arguments{
  \item{object}{an object of class PvInd}

  \item{cov}{a character vector indicating which covariate
  have to be used in the analyses. As current
  implementation of the lasso (glmnet) does not Such
  variables are first use in a logistic regression model
  and the residual is then used as the outcome in the
  penalized regression}

  \item{aeLab}{The label of the adverse event to be
  regressed. By default, the \code{pvPen} will regress all
  ae, one at a time and this can be very time consuming and
  RAM demanding.}

  \item{detectCriter}{Can be either BIC or AIC. These
  criteria are used to select the final model.}

  \item{posConst}{If TRUE (default), the regression
  coefficients are constrained to be non non negative, this
  to ensure that the final model only contains drugs
  "increasing" the risk of a given ae.}

  \item{nDrugMax}{Maximum number of drugs to be included in
  the model. In addition to be computationnaly intensive,
  odels with to many drugs are likely to be very unstable.}

  \item{parallel}{Whether parallel should be used in order
  to speed up the calculation. Be careful as it will be
  more RAM demanding. Only available for linux or Mac Os}

  \item{nCore}{If parallel is true, allows to specify the
  number of cores to be used for parallel computing.}

  \item{...}{Further arguments to be passed to other
  functions (None for the moment)}
}
\description{
  pvPen implements a detection strategy based on the use of
  the Lasso regression along with BIC type detection
  criteria. The methodology follows two main steps:
  \itemize{ \item Lasso regression (with a positive
  constraint on the regression coefficients) for a grid of
  constraint parameters \item selection on the drug
  associated with the ae of interest according to the BIC.
  } This function heavily relies on the highly efficient
  \code{glmnet} package. \code{pvPen} is computationnaly
  highly demanding. It is also strongly advised to use a
  linux or a mac computer in order to use several cores.
  Note also that this regression approach should be used
  with ae and drugs both having a reasonable number of
  reports.
}

