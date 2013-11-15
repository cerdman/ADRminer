\docType{class}
\encoding{UTF-8}
\name{pvInd-class}
\alias{getAe}
\alias{getAeMargin}
\alias{getAeMargin,pvInd-method}
\alias{getAe,pvInd-method}
\alias{getCov}
\alias{getCov,pvInd-method}
\alias{getDrug}
\alias{getDrugMargin}
\alias{getDrugMargin,pvInd-method}
\alias{getDrug,pvInd-method}
\alias{names,pvInd-method}
\alias{pvInd-class}
\alias{$<-,pvInd-method}
\alias{$,pvInd-method}
\alias{shows,pvInd-method}
\title{pvInd class}
\description{
  pvInd class
}
\section{Methods}{
  \describe{ \item{names}{\code{signature(x = "pvInd")}:
  returns the names of the slots of the object.}
  \item{show}{\code{signature(x = "pvInd")}: printing of
  the object.} \item{$}{\code{signature(x = "pvInd")}:
  similar to the @ operator; used to access the content of
  slots of the object.} \item{$<-}{\code{signature(x =
  "pvInd")}: similar to the @ operator; used to replace the
  content of slots of the object.}
  \item{getDrug}{\code{signature(x = "pvInd")}: returns the
  names of the slots of the object.}
  \item{getAe}{\code{signature(x = "pvInd")}: returns the
  adverse event matrix.}
  \item{getDrugMargin}{\code{signature(x = "pvInd")}:
  returns the Drug margin counts.}
  \item{getAeMargin}{\code{signature(x = "pvInd")}: returns
  the adverse event margin counts.}
  \item{getCov}{\code{signature(x = "pvInd")}: returns the
  covariate data.frame.} }
}
\author{
  Youness Ergaibi & Ismaïl Ahmed
}

