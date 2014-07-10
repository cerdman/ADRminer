\encoding{UTF-8}
\name{prr}
\alias{prr}
\title{Proportional Reporting Ratio}
\usage{
  prr(object, rr0 = 1, assocMeasure = c("pvalue", "Lb95"),
    detectCriter = c("FDR", "nbSig", "assocMeasure"),
    criterThres = 0.05, nMin = 1, strat = NULL, allRes = F)
}
\arguments{
  \item{object}{An object of class PvInd or PvCont.}

  \item{rr0}{The tested relative risk. By default, rr0=1.}

  \item{assocMeasure}{Statistic used to order the
  drug-events pairs: pvalue: P-value Lb95: Lower bound of
  the 95\% two sided confidence interval of log(prr).}

  \item{detectCriter}{Decision rule for the signal
  generation based on: \itemize{ \item FDR: a prespecified
  value (\code{detectCriter}) for the false discovery rate
  (Default value) \item nbSig: a prespecified number
  (\code{detectCriter}) of signals \item assocMeasure: a
  prespecified value (\code{detectCriter}) for the Ranking
  statistic (\code{assocMeasure}) }}

  \item{criterThres}{Threshold for \code{detectCriter}. Ex
  0.05 for FDR.}

  \item{nMin}{Minimum number of spontaneous reports for a
  drug-event pair to be potentially considered as a signal.
  By default, \code{nMin}=1.}

  \item{strat}{Character vector containing the name of the
  covariates to be used for a stratified analysis}

  \item{allRes}{Logical value indicating whether all
  drug-event combination results must be provided.}
}
\value{
  allSig Data.frame summarizing the results for all
  drug-event combinations with at least \code{nMin}
  spontaneous reports ordered by \code{assocMeasure}. This
  output is provided if \code{allRes=TRUE}. Operating
  characteristics are estimated according to Ahmed et al
  (Biometrics 2010).

  sig Same as \code{allSig} but restricted to the list of
  generated signals.

  nSig Number of generated signals.

  call Arguments entered when calling \code{prr}.
}
\description{
  Proportional Reporting Ratio initially proposed by Evans
  et al. (2001) extended to the multiple comparison
  framework. Note that the computed variance is different
  from the one used in van Puijenbroek et al. (2002). The
  rule proposed by Evans et al. is not implemeted.
}
\author{
  Youness Ergaibi & Ismaïl Ahmed
}
\keyword{prr}

