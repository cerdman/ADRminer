\docType{methods}
\encoding{UTF-8}
\name{pvInd}
\alias{names,pvInd}
\alias{pvInd}
\alias{pvInd,data.frame,data.frame-method}
\alias{pvInd,data.frame,missing-method}
\alias{shows,pvInd}
\title{pvInd}
\arguments{
  \item{adr}{a data.frame with three columns:}

  \item{cov}{data.frame which contains individual
  covariates that can be used for stratified analysis. It
  is assumed that the first column of the data.frame
  corresponds to the spontaneous report identifier.}

  \item{dMarginLim}{Minimum number of spontaneous reports
  for a drug to be included (default is 1)}

  \item{aeMarginLim}{Minimum number of spontaneous reports
  for an adverse event to be included (default is 1)}
}
\value{
  a pvInd object
}
\description{
  \code{pvInd} is used to convert raw data (\code{adr} and
  ) \code{cov}) into an pvInd object that can be used in
  the signal detection method functions.
}
\author{
  Ismaïl Ahmed
}
\keyword{pvInd}

