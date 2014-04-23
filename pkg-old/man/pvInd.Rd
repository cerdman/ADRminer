\docType{methods}
\encoding{UTF-8}
\name{pvInd}
\alias{pvInd}
\alias{pvInd,data.frame,data.frame-method}
\alias{pvInd,data.frame,missing-method}
\title{pvInd}
\arguments{
  \item{adr}{a data.frame with three columns:}

  \item{cov}{data.frame which contains individual
  covariates that can be used for stratified analysis. It
  is assumed that the first column of the data.frame
  corresponds to the spontaneous report identifier.}

  \item{drugMarginLim}{Minimum number of spontaneous
  reports for a drug to be included (default is 1)}

  \item{aeMarginLim}{Minimum number of spontaneous reports
  for an adverse event to be included (default is 1)}
}
\value{
  a pvInd object
}
\description{
  This function (which takes the same name that its
  corresponding formal class \code{\link{pvInd-class}})
  should be used to convert a data.frame containing
  individual spontaneous reports (and optionnaly
  supplementary individual information) into an pvInd
  object.
}
\author{
  Ismaïl Ahmed
}
\seealso{
  \code{\link{pvInd-class}}
}
\keyword{pvInd}

