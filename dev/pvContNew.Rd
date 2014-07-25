\encoding{UTF-8}
\name{pvContNew}
\alias{pvContNew}
\title{pvContNew}
\arguments{
  \item{object}{a PvInd object or a data.frame. The latter
  should have 3 columns, i.e. the drug labels, the ae
  labels and the corresponding number of spontaneous
  reports.}

  \item{strat}{Character vector containing the name of the
  covariates to be used for a stratified analysis.
  \code{strat} is only meaningful when used along with a
  PvInd object.}
}
\value{
  A list of pvCont objects of class
  \code{\link{PvContStrat}}.
}
\description{
  pvContNew converts a pvInd object into PvCont object. If
  a \code{strat} is provided, it outputs a list of pvCont
  objects (pvContStrat object) for each combination of
  \code{strat}
}
\author{
  Ismaïl Ahmed
}
\keyword{pvContNew}

