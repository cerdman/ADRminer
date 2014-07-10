\name{bcpnn}
\alias{bcpnn}
\title{Bayesian confidence propagation neural network}
\usage{
  bcpnn(object, rr0 = 1, assocMeasure = c("post.H0", "LB"),
    detectCriter = c("FDR", "nbSig", "assocMeasure"),
    criterThres = 0.05, nMin = 1, MC = FALSE,
    NB.MC = 10000, allRes = TRUE)
}
\arguments{
  \item{object}{is an object(pvInd or pvCont)}

  \item{rr0}{is the Value of the tested risk. By default,
  RR0=1.}

  \item{assocMeasure}{is a Statistic used for ranking the
  couples: post.H0:Posterior probability of the null
  hypothesis LB:2.5% quantile of the posterior distribution
  of IC.}

  \item{detectCriter}{is a Decision rule for the signal
  generation based on: FDR: (Default value) nbSig:Number of
  signals assocMeasure:Ranking statistic. See
  \code{assocMeasure}}

  \item{criterThres}{is a Threshold for
  \code{detectCriter}. Ex 0.05 for FDR.}

  \item{nMin}{is a Minimum number of notifications for a
  couple to be potentially considered as a signal. By
  default, \code{nMin=1}.}

  \item{NB.MC}{If MC=TRUE, NB.MC indicates the number of
  Monte Carlo simulations to be done}

  \item{MC}{If MC=TRUE, the statistic of interest (see
  RANKSTAT) is calculated by Monte Carlo simulations which
  can be very long. If MC=FALSE, IC is approximated by a
  normal distribution (which can be very crude for small
  counts).}

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
  Bayesian confidence propagation neural network (Bate et
  al. 1998, Noren et al. 2006) extended to the multiple
  comparison framework.
}
\details{
  \code{bcpnn} function
}
\author{
  Youness Ergaibi & Ismaïl Ahmed
}
\keyword{bcpnn}

