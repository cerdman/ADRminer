\name{rfet}
\alias{rfet}
\title{Reporting Fisher's Exact Test}
\usage{
  rfet(object, or0 = 1,
    detectCriter = c("FDR", "nbSig", "pvalue"),
    criterThres = 0.05, MID.PVAL = FALSE, nMin = 1,
    allRes = F)
}
\arguments{
  \item{object}{is an object(pvInd or pvCont)}

  \item{or0}{is a value of the tested odds ratio. By
  default or0=1}

  \item{detectCriter}{is a Decision rule for the signal
  generation based on: FDR: (Default value) nbSig:Number of
  signals pvalue:Ranking statistic. See \code{pvalue}}

  \item{criterThres}{is a Threshold for
  \code{detectCriter}. Ex 0.05 for FDR.}

  \item{MID.PVAL}{the statistic of interest becomes the
  mid-P-values instead of the P-values resulting from the
  Fisher's exact test. By default \code{MID.PVAL=FALSE}.}

  \item{nMin}{is a Minimum number of notifications for a
  couple to be potentially considered as a signal. By
  default, \code{nMin=1}.}

  \item{allRes}{is a logical value}
}
\value{
  allSig Data.frame summarizing the results of all couples
  with at least \code{nMin} notifications ordered by
  \code{assocMeasure}. It contains notably the labels, the
  cell counts, the expected counts n1. * n.1 / N

  sig Same Data.frame as \code{allSig} but restricted to
  the list of generated signals.

  nbSig Number of generated signals.

  call Parameters entered in the function.
}
\description{
  Reporting Odds Ratio proposed by van Puijenbroak et al.
  (2002) extended to the multiple comparison framework.
}
\details{
  \code{rfet} function
}
\author{
  Youness Ergaibi & Ismaïl Ahmed
}
\keyword{RFET}

