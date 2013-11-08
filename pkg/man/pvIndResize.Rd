\encoding{UTF-8}
\name{pvIndResize}
\alias{pvIndResize}
\title{pvIndResize}
\usage{
  pvIndResize(pvIndObj, aeMarginMin = 1, dMarginMin = 1)
}
\arguments{
  \item{pvIndObj}{an object of class "pvInd" (see
  \code{\link{pvInd-class}} and \code{\link{pvInd}})}

  \item{aeMarginMin}{Minimum number of reports in which an
  ae must be involved in order to be kept}

  \item{dMarginMin}{Minimum number of reports in which a
  drug must be involved in order to be kept}
}
\value{
  an object of class pvInd
}
\description{
  \code{pvIndResize} makes it possible to remove drug(s)
  and/or ae(s) whose margin are less than a given number.
}
\author{
  Ismaïl Ahmed
}

