\encoding{UTF-8}
\name{cv.glmnet.pvInd}
\alias{cv.glmnet.pvInd}
\title{Cross-validation for glmnet and pvInd object}
\usage{
  cv.glmnet.pvInd(object, aeId = "all", covId = NULL,
    posConst = TRUE, nDrugMax = 20, parallel = FALSE, ...)
}
\arguments{
  \item{object}{an object of class pvInd}

  \item{aeId}{The label of the adverse event to be
  regressed. By default, the \code{cv.glmnet} will regress
  all ae, one at a time and this can be very time consuming
  and RAM demanding. This parameter can be either the name
  of the ae(s), or the index of the ae column}

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

  \item{...}{Further arguments to be passed to the
  cv.glmnet function (not in use for now)}
}
\description{
  This function is mainly a wrapper of the cv.glmnet
  function (glmnet package) for pvInd object
}

